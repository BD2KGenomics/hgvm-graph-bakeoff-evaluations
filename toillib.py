"""
toillib.py: useful extras for Toil scripts.

Includes real-time-logging, input retrieval from file/S3/Azure, and output
deposit to file/S3/Azure

The only problem is you can't use it since it won't be installed on your nodes
"""


import sys, os, os.path, json, collections, logging, logging.handlers
import SocketServer, struct, socket, threading, tarfile, shutil

# We need some stuff in order to have Azure
try:
    import azure
    # Make sure to get the 0.11 BlobService, in case the new azure storage
    # module is also installed.
    from azure.storage import BlobService
    import toil.jobStores.azureJobStore
    have_azure = True
except ImportError:
    have_azure = False
    pass
    

def robust_makedirs(directory):
    """
    Make a directory when other nodes may be trying to do the same on a shared
    filesystem.
    
    """

    if not os.path.exists(directory):
        try:
            # Make it if it doesn't exist
            os.makedirs(directory)
        except OSError:
            # If you can't make it, maybe someone else did?
            pass
            
    # Make sure it exists and is a directory
    assert(os.path.exists(directory) and os.path.isdir(directory))

class LoggingDatagramHandler(SocketServer.DatagramRequestHandler):
    """
    Receive logging messages from the jobs and display them on the master.
    
    Uses length-prefixed JSON message encoding.
    """
    
    def handle(self):
        """
        Handle messages coming in over self.connection.
        
        Messages are 4-byte-length-prefixed JSON-encoded logging module records.
        """
        
        while True:
            # Loop until we run out of messages
        
            # Parse the length
            length_data = self.rfile.read(4)
            if len(length_data) < 4:
                # The connection was closed, or we didn't get enough data
                # TODO: complain?
                break
                
            # Actually parse the length
            length = struct.unpack(">L", length_data)[0]
            
            # This is where we'll put the received message
            message_parts = []
            length_received = 0
            while length_received < length:
                # Keep trying to get enough data
                part = self.rfile.read(length - length_received)
                
                length_received += len(part)
                message_parts.append(part)
                
            # Stitch it all together
            message = "".join(message_parts)

            try:
            
                # Parse it as JSON
                message_attrs = json.loads(message)
                
                # Fluff it up into a proper logging record
                record = logging.makeLogRecord(message_attrs)
            except:
                logging.error("Malformed record")
                
            # TODO: do log level filtering
            logging.getLogger("remote").handle(record)
            
class JSONDatagramHandler(logging.handlers.DatagramHandler):
    """
    Send logging records over UDP serialized as JSON.
    """
    
    def makePickle(self, record):
        """
        Actually, encode the record as length-prefixed JSON instead.
        """
        
        json_string = json.dumps(record.__dict__)
        length = struct.pack(">L", len(json_string))
        
        return length + json_string
        
class RealTimeLogger(object):
    """
    All-static class for getting a logger that logs over UDP to the master.
    """
    
    # Also the logger
    logger = None
    
    # The master keeps a server and thread
    logging_server = None
    server_thread = None
  
    @classmethod
    def start_master(cls):
        """
        Start up the master server and put its details into the options
        namespace.
        
        """
        
        logging.basicConfig(level=logging.INFO)
    
        # Start up the logging server
        cls.logging_server = SocketServer.ThreadingUDPServer(("0.0.0.0", 0),
            LoggingDatagramHandler)
            
        # Set up a thread to do all the serving in the background and exit when we
        # do
        cls.server_thread = threading.Thread(
            target=cls.logging_server.serve_forever)
        cls.server_thread.daemon = True
        cls.server_thread.start()
        
        # Set options for logging in the class and the options namespace
        # Save them in the environment so they get sent out to jobs
        os.environ["RT_LOGGING_HOST"] = socket.getfqdn()
        os.environ["RT_LOGGING_PORT"] = str(
            cls.logging_server.server_address[1])
        
        
    @classmethod
    def stop_master(cls):
        """
        Stop the server on the master.
        
        """
        
        cls.logging_server.shutdown()
        cls.server_thread.join()
  
    @classmethod
    def get(cls):
        """
        Get the logger that logs to master.
        
        Note that if the master logs here, you will see the message twice,
        since it still goes to the normal log handlers too.
        """
        
        if cls.logger is None:
            # Only do the setup once, so we don't add a handler every time we
            # log
            cls.logger = logging.getLogger('realtime')
            cls.logger.setLevel(logging.INFO)
            
            if (os.environ.has_key("RT_LOGGING_HOST") and
                os.environ.has_key("RT_LOGGING_PORT")):
                # We know where to send messages to, so send them.
            
                cls.logger.addHandler(JSONDatagramHandler(
                    os.environ["RT_LOGGING_HOST"],
                    int(os.environ["RT_LOGGING_PORT"])))
        
        return cls.logger

def write_global_directory(file_store, path, cleanup=False):
    """
    Write the given directory into the file store, and return an ID that can be
    used to retrieve it. Writes the files in the directory and subdirectories
    into a tar file in the file store.

    Does not preserve the name or permissions of the given directory (only of
    its contents).

    If cleanup is true, directory will be deleted from the file store when this
    job and its follow-ons finish.
    
    """
    
    with file_store.writeGlobalFileStream(cleanup=cleanup) as (file_handle,
        file_id):
        # We have a stream, so start taring into it
    
        with tarfile.open(fileobj=file_handle, mode="w|gz") as tar:
            # Open it for streaming-only write (no seeking)
            
            # We can't just add the root directory, since then we wouldn't be
            # able to extract it later with an arbitrary name.
            
            for file_name in os.listdir(path):
                # Add each file in the directory to the tar, with a relative
                # path
                tar.add(os.path.join(path, file_name), arcname=file_name)
                
        # Spit back the ID to use to retrieve it
        return file_id
        
def read_global_directory(file_store, directory_id, path):
    """
    Reads a directory with the given tar file id from the global file store and
    recreates it at the given path.
    
    The given path, if it exists, must be a directory.
    
    Do not use to extract untrusted directories, since they could sneakily plant
    files anywhere on the filesystem.
    
    """
    
    # Make the path
    robust_makedirs(path)
    
    with file_store.readGlobalFileStream(directory_id) as file_handle:
        # We need to pull files out of this tar stream
    
        with tarfile.open(fileobj=file_handle, mode="r|*") as tar:
            # Open it for streaming-only read (no seeking)
            
            # We need to extract the whole thing into that new directory
            tar.extractall(path)
            

class IOStore(object):
    """
    A class that lets you get your input files and save your output files
    to/from a local filesystem, Amazon S3, or Microsoft Azure storage
    transparently.
    
    This is the abstract base class; other classes inherit from this and fill in
    the methods.
    
    """
    
    def __init__(self):
        """
        Make a new IOStore
        """
        
        raise NotImplementedError()
            
    def read_input_file(self, input_path, local_path):
        """
        Read an input file from wherever the input comes from and send it to the
        given path.
        
        """
        
        raise NotImplementedError()
        
    def list_input_directory(self, input_path, recursive=False):
        """
        Yields each of the subdirectories and files in the given input path.
        
        If recursive is false, yields files and directories in the given
        directory. If recursive is true, yields all files contained within the
        current directory, recursively, but does not yield folders.
        
        
        Gives relative file/directory names.
        
        """
        
        raise NotImplementedError()
    
    def write_output_file(self, local_path, output_path):
        """
        Save the given local file to the given output path. No output directory
        needs to exist already.
        
        """
        
        raise NotImplementedError()
        
    def exists(self, path):
        """
        Returns true if the given input or output file exists in the store
        already.
        
        """
        
        raise NotImplementedError()
        
    @staticmethod
    def get(store_string):
        """
        Get a concrete IOStore created from the given connection string.
        
        Valid formats are just like for a Toil JobStore, except with container
        names being specified on Azure.
        
        Formats:
        
        /absolute/filesystem/path
        
        ./relative/filesystem/path
        
        file:filesystem/path
        
        aws:region:bucket (TODO)
        
        aws:region:bucket/path/prefix (TODO)
        
        azure:account:container (instead of a container prefix) (gets keys like
        Toil)
        
        azure:account:container/path/prefix (trailing slash added automatically)
        
        """
        
        # Code adapted from toil's common.py loadJobStore()
        
        if store_string[0] in "/.":
            # Prepend file: tot he path
            store_string = "file:" + store_string

        try:
            # Break off the first colon-separated piece.
            store_type, store_arguments = store_string.split(":", 1)
        except ValueError:
            # They probably forgot the . or /
            raise RuntimeError("Incorrect IO store specification {}. "
                "Local paths must start with . or /".format(store_string))

        if store_type == "file":
            return FileIOStore(store_arguments)
        elif store_type == "aws":
            # Break out the AWS arguments
            region, name_prefix = store_arguments.split(":", 1)
            raise NotImplementedError("AWS IO store not written")
        elif store_type == "azure":
            # Break out the Azure arguments. TODO: prefix in the container.
            account, container = store_arguments.split(":", 1)
            
            if "/" in container:
                # Split the container from the path
                container, path_prefix = container.split("/", 1)
            else:
                # No path prefix
                path_prefix = ""
            
            return AzureIOStore(account, container, path_prefix)
        else:
            raise RuntimeError("Unknown IOStore implementation {}".format(
                store_type))
        
        
            
class FileIOStore(IOStore):
    """
    A class that lets you get input from and send output to filesystem files.
    
    """
    
    def __init__(self, path_prefix=""):
        """
        Make a new FileIOStore that just treats everything as local paths,
        relative to the given prefix.
        
        """
        
        self.path_prefix = path_prefix
        
    def read_input_file(self, input_path, local_path):
        """
        Get input from the filesystem.
        """
        
        RealTimeLogger.get().debug("Loading {} from FileIOStore in {}".format(
            input_path, self.path_prefix))
        
        # Make a symlink to grab things
        os.symlink(os.path.abspath(os.path.join(self.path_prefix, input_path)),
            local_path)
        
    def list_input_directory(self, input_path, recursive=False):
        """
        Loop over directories on the filesystem.
        """
        
        RealTimeLogger.get().info("Enumerating {} from "
            "FileIOStore in {}".format(input_path, self.path_prefix))
        
        for item in os.listdir(os.path.join(self.path_prefix, input_path)):
            if(recursive and os.isdir(item)):
                # Recurse on this
                for subitem in self.list_input_directory(
                    os.path.join(input_path, item), recursive):
                    
                    # Make relative paths include this directory anme and yield
                    # them
                    yield os.path.join(item, subitem)
            else:
                # This isn't a directory or we aren't being recursive
                yield item
    
    def write_output_file(self, local_path, output_path):
        """
        Write output to the filesystem
        """

        RealTimeLogger.get().debug("Saving {} to FileIOStore in {}".format(
            output_path, self.path_prefix))

        # What's the real outptu path to write to?
        real_output_path = os.path.join(self.path_prefix, output_path)

        # What directory should this go in?
        parent_dir = os.path.split(real_output_path)[0]
            
        if parent_dir != "":
            # Make sure the directory it goes in exists.
            robust_makedirs(parent_dir)
        
        # These are small so we just make copies
        shutil.copy2(local_path, real_output_path)
        
    def exists(self, path):
        """
        Returns true if the given input or output file exists in the file system
        already.
        
        """
        
        return os.path.exists(os.path.join(self.path_prefix, path))
            
class AzureIOStore(IOStore):
    """
    A class that lets you get input from and send output to Azure Storage.
    
    """
    
    def __init__(self, account_name, container_name, name_prefix=""):
        """
        Make a new AzureIOStore that reads from and writes to the given
        container in the given account, adding the given prefix to keys. All
        paths will be interpreted as keys or key prefixes.
        
        If the name prefix does not end with a trailing slash, and is not empty,
        one will be added automatically.
        
        Account keys are retrieved from the AZURE_ACCOUNT_KEY environment
        variable or from the ~/.toilAzureCredentials file, as in Toil itself.
        
        """
        
        # Make sure azure libraries actually loaded
        assert(have_azure)
        
        self.account_name = account_name
        self.container_name = container_name
        self.name_prefix = name_prefix
        
        if self.name_prefix != "" and not self.name_prefix.endswith("/"):
            # Make sure it has the trailing slash required.
            self.name_prefix += "/"
        
        # Sneak into Toil and use the same keys it uses
        self.account_key = toil.jobStores.azureJobStore._fetchAzureAccountKey(
            self.account_name)
            
        # This will hold out Azure blob store connection
        self.connection = None
        
    def __getstate__(self):
        """
        Return the state to use for pickling. We don't want to try and pickle
        an open Azure connection.
        """
     
        return (self.account_name, self.account_key, self.container_name, 
            self.name_prefix)
        
    def __setstate__(self, state):
        """
        Set up after unpickling.
        """
        
        self.account_name = state[0]
        self.account_key = state[1]
        self.container_name = state[2]
        self.name_prefix = state[3]
        
        self.connection = None
        
    def __connect(self):
        """
        Make sure we have an Azure connection, and set one up if we don't.
        """
        
        if self.connection is None:
            RealTimeLogger.get().debug("Connecting to account {}, using "
                "container {} and prefix {}".format(self.account_name,
                self.container_name, self.name_prefix))
        
            # Connect to the blob service where we keep everything
            self.connection = BlobService(
                account_name=self.account_name, account_key=self.account_key)
            
            
    def read_input_file(self, input_path, local_path):
        """
        Get input from Azure.
        """
        
        self.__connect()
        
        
        RealTimeLogger.get().debug("Loading {} from AzureIOStore".format(
            input_path))
        
        # Download the blob. This is known to be synchronous, although it can
        # call a callback during the process.
        self.connection.get_blob_to_path(self.container_name,
            self.name_prefix + input_path, local_path)
            
    def list_input_directory(self, input_path, recursive=False):
        """
        Loop over fake /-delimited directories on Azure. The prefix may or may
        not not have a trailing slash; if not, one will be added automatically.
        
        Returns the names of files and fake directories in the given input fake
        directory, non-recursively.
        
        """
        
        self.__connect()
        
        RealTimeLogger.get().info("Enumerating {} from AzureIOStore".format(
            input_path))
        
        # Work out what the directory name to list is
        fake_directory = self.name_prefix + input_path
        
        if fake_directory != "" and not fake_directory.endswith("/"):
            # We have a nonempty prefix, and we need to end it with a slash
            fake_directory += "/"
        
        # This will hold the marker that we need to send back to get the next
        # page, if there is one. See <http://stackoverflow.com/a/24303682>
        marker = None
        
        # This holds the subdirectories we found; we yield each exactly once if
        # we aren't recursing.
        subdirectories = set()
        
        while True:
        
            # Get the results from Azure. We skip the delimiter since it doesn't
            # seem to have the placeholder entries it's suppsoed to.
            result = self.connection.list_blobs(self.container_name, 
                prefix=fake_directory, marker=marker)
                
            for blob in result:
                # Yield each result's blob name, but directory names only once
                
                # Drop the common prefix
                relative_path = blob.name[len(fake_directory):]
                
                if (not recursive) and "/" in relative_path:
                    # We found a file in a subdirectory, and we aren't supposed
                    # to be recursing.
                    subdirectory, _ = relative_path.split("/", 1)
                    
                    if subdirectory not in subdirectories:
                        # It's a new subdirectory. Yield and remember it
                        subdirectories.add(subdirectory)
                        
                        yield subdirectory
                else:
                    # We found an actual file  
                    yield relative_path
                
            # Save the marker
            marker = result.next_marker
                
            if not marker:
                break 
    
    def write_output_file(self, local_path, output_path):
        """
        Write output to Azure. Will create the container if necessary.
        """
        
        self.__connect()
        
        RealTimeLogger.get().debug("Saving {} to AzureIOStore".format(
            output_path))
        
        try:
            # Make the container
            self.connection.create_container(self.container_name)
        except azure.WindowsAzureConflictError:
            # The container probably already exists
            pass
        
        # Upload the blob (synchronously)
        # TODO: catch no container error here, make the container, and retry
        self.connection.put_block_blob_from_path(self.container_name,
            self.name_prefix + output_path, local_path)
            
    def exists(self, path):
        """
        Returns true if the given input or output file exists in Azure already.
        
        """
        
        self.__connect()
        
        marker = None
        
        while True:
        
            # Get the results from Azure.
            result = self.connection.list_blobs(self.container_name, 
                prefix=self.name_prefix + path, marker=marker)
                
            for blob in result:
                # Look at each blob
                
                if blob.name == self.name_prefix + path:
                    # Found it
                    return True
                
            # Save the marker
            marker = result.next_marker
                
            if not marker:
                break 
        
        return False

    





















