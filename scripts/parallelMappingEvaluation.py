#!/usr/bin/env python2.7
"""
parallelMappingEvaluation.py: Run the mapping evaluation on all the servers in
parallel.

BAM files with reads must have been already downloaded.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading

from toil.job import Job

###BEGIN TOILLIB

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
        
        logging.basicConfig(level=logging.DEBUG)
    
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
            cls.logger.setLevel(logging.DEBUG)
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
        
    def list_input_directory(self, input_path):
        """
        Yields each of the subdirectories and files in the given input path, non-
        recursively.
        
        Gives bare file/directory names with no paths.
        
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
        
        RealTimeLogger.get().info("Loading {} from FileIOStore in {}".format(
            input_path, self.path_prefix))
        
        # Make a symlink to grab things
        os.symlink(os.path.abspath(os.path.join(self.path_prefix, input_path)),
            local_path)
        
    def list_input_directory(self, input_path):
        """
        Loop over directories on the filesystem.
        """
        
        RealTimeLogger.get().info("Enumerating {} from "
            "FileIOStore in {}".format(input_path, self.path_prefix))
        
        for item in os.listdir(os.path.join(self.path_prefix, input_path)):
            yield item
    
    def write_output_file(self, local_path, output_path):
        """
        Write output to the filesystem
        """

        RealTimeLogger.get().info("Saving {} to FileIOStore in {}".format(
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
            RealTimeLogger.get().info("Connecting to account {}, using "
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
        
        
        RealTimeLogger.get().info("Loading {} from AzureIOStore".format(
            input_path))
        
        # Download the blob. This is known to be synchronous, although it can
        # call a callback during the process.
        self.connection.get_blob_to_path(self.container_name,
            self.name_prefix + input_path, local_path)
            
    def list_input_directory(self, input_path):
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
        
        # This holds the subdirectories we found; we yield each exactly once.
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
                
                if "/" in relative_path:
                    # We found a file in a subdirectory
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
        
        RealTimeLogger.get().info("Saving {} to AzureIOStore".format(
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

###END TOILLIB

def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("server_list", type=argparse.FileType("r"),
        help="TSV file continaing <region>\t<url> lines for servers to test")
    parser.add_argument("sample_store",
        help="sample input IOStore with <region>/<sample>/<sample>.bam.fq")
    parser.add_argument("out_store",
        help="output IOStore to create and fill with alignments and stats")
    parser.add_argument("--server_version", default="v0.6.g",
        help="server version to add to URLs")
    parser.add_argument("--sample_limit", type=int, default=1, 
        help="number of samples to use")
    parser.add_argument("--edge_max", type=int, default=0, 
        help="maximum edges to cross in index")
    parser.add_argument("--kmer_size", type=int, default=10, 
        help="size of kmers to use in indexing and mapping")
    parser.add_argument("--bin_url",
        default="https://hgvm.blob.core.windows.net/hgvm-bin",
        help="URL to download sg2vg and vg binaries from, without Docker")
    parser.add_argument("--overwrite", default=False, action="store_true",
        help="overwrite existing result files")
    parser.add_argument("--reindex", default=False, action="store_true",
        help="don't re-use existing indexed graphs")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    

def run_all_alignments(job, options):
    """
    For each server listed in the server_list tsv, kick off child jobs to
    align and evaluate it.

    """
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    sample_store = IOStore.get(options.sample_store)
    out_store = IOStore.get(options.out_store)
    
    # Retrieve binaries we need
    RealTimeLogger.get().info("Retrieving binaries from {}".format(
        options.bin_url))
    bin_dir = "{}/bin".format(job.fileStore.getLocalTempDir())
    robust_makedirs(bin_dir)
    subprocess.check_call(["wget", "{}/sg2vg".format(options.bin_url),
        "-O", "{}/sg2vg".format(bin_dir)])
    subprocess.check_call(["wget", "{}/vg".format(options.bin_url),
        "-O", "{}/vg".format(bin_dir)])
        
    # Make them executable
    os.chmod("{}/sg2vg".format(bin_dir), 0o744)
    os.chmod("{}/vg".format(bin_dir), 0o744)
    
    # Upload the bin directory to the file store
    bin_dir_id = write_global_directory(job.fileStore, bin_dir,
        cleanup=True)
    
    # Make sure we skip the header
    is_first = True
    
    for line in options.server_list:
        if is_first:
            # This is the header, skip it.
            is_first = False
            continue
        
        # We need to read each non-header line
        
        # Break it into its fields
        parts = line.split("\t")
        
        if parts[0].startswith("#"):
            # Skip comments
            continue
            
        if parts[0].startswith("\n"):
            # Skip newlines
            continue
            
        # Pull out the first 3 fields
        region, url, generator = parts[0:3]
        
        # We cleverly just split the lines out to different nodes
        job.addChildJobFn(run_region_alignments, options, bin_dir_id, region,
            url, cores=16, memory="100G", disk="50G")
            
        # Say what we did
        RealTimeLogger.get().info("Running child for {}".format(parts[1]))
        

def run_region_alignments(job, options, bin_dir_id, region, url):
    """
    For the given region, download, index, and then align to the given graph.
    
    """
    
    RealTimeLogger.get().info("Running on {} for {}".format(url, region))
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    sample_store = IOStore.get(options.sample_store)
    out_store = IOStore.get(options.out_store)
    
    # Download the binaries
    bin_dir = "{}/bin".format(job.fileStore.getLocalTempDir())
    read_global_directory(job.fileStore, bin_dir_id, bin_dir)
    
    # Get graph basename (last URL component) from URL
    basename = re.match(".*/(.*)/$", url).group(1)
        
    # Get graph name (without region and its associated dash) from basename
    graph_name = basename.replace("-{}".format(region), "").replace(
        "{}-".format(region), "")
    
    # Where do we look for samples for this region in the input?
    region_dir = region.upper()
    
    # What samples do we do? List input sample names up to the given limit.
    input_samples = list(sample_store.list_input_directory(
        region_dir))[:options.sample_limit]
    
    # Work out the directory for the alignments to be dumped in in the output
    alignment_dir = "alignments/{}/{}".format(region, graph_name)
    
    # Also for statistics
    stats_dir = "stats/{}/{}".format(region, graph_name)
    
    # What samples haven't been done yet and need doing
    samples_to_run = []
    
    for sample in input_samples:
        # Split out over each sample
        
        # What's the file that has to exist for us to not re-run it?
        stats_file_key = "{}/{}.json".format(stats_dir, sample)
        
        if (not options.overwrite) and out_store.exists(stats_file_key):
            # This is already done.
            RealTimeLogger.get().info("Skipping completed alignment of "
                "{} to {} {}".format(sample, graph_name, region))
            continue
        else:
            # We need to run this sample
            samples_to_run.append(sample)
            
    if len(samples_to_run) == 0 and not options.reindex:
        # Don't bother indexing the graph if all the samples are done, and we
        # didn't explicitly ask to do it.
        RealTimeLogger.get().info("Nothing to align to {}".format(basename))
        return
    
    # Make the real URL with the version
    versioned_url = url + options.server_version
    
    # Where will the indexed graph go in the output
    index_key = "indexes/{}/{}.tar.gz".format(region, graph_name)
    
    if (not options.reindex) and out_store.exists(index_key):
        # See if we have an index already available in the output store from a
        # previous run
        
        RealTimeLogger.get().info("Retrieving indexed {} graph from output "
            "store".format(basename))
            
        # Download the pre-made index directory
        tgz_file = "{}/index.tar.gz".format(job.fileStore.getLocalTempDir())
        out_store.read_input_file(index_key, tgz_file)
        
        # Save it to the global file store and keep around the ID.
        # Will be compatible with read_global_directory
        index_dir_id = job.fileStore.writeGlobalFile(tgz_file, cleanup=True)
        
    else:
        # Download the graph, build the index, and store it in the output store
    
        # Work out where the graph goes
        # it will be graph.vg in here
        graph_dir = "{}/graph".format(job.fileStore.getLocalTempDir())
        robust_makedirs(graph_dir)
        
        graph_filename = "{}/graph.vg".format(graph_dir)
        
        # Download and fix up the graph with this ugly subprocess pipeline
        # sg2vg "${URL}" -u | vg view -Jv - | vg mod -X 100 - | 
        # vg ids -s - > "graphs/${BASENAME}.vg"
        
        with open(graph_filename, "w") as output_file:
        
            RealTimeLogger.get().info("Downloading {} to {}".format(
                versioned_url, graph_filename))
        
            # Hold all the popen objects we need for this
            tasks = []
            
            # Do the download
            tasks.append(subprocess.Popen(["{}/sg2vg".format(bin_dir),
                versioned_url, "-u"], stdout=subprocess.PIPE))
            
            # Pipe through zcat
            tasks.append(subprocess.Popen(["{}/vg".format(bin_dir), "view",
                "-Jv", "-"], stdin=tasks[-1].stdout, stdout=subprocess.PIPE))
            
            # And cut
            tasks.append(subprocess.Popen(["{}/vg".format(bin_dir), "mod",
                "-X100", "-"], stdin=tasks[-1].stdout, stdout=subprocess.PIPE))
                
            # And uniq
            tasks.append(subprocess.Popen(["{}/vg".format(bin_dir), "ids", "-s",
                "-"], stdin=tasks[-1].stdout, stdout=output_file))
                
            # Did we make it through all the tasks OK?
            for task in tasks:
                if task.wait() != 0:
                    raise RuntimeError("Pipeline step returned {}".format(
                        task.returncode))
        
        # Now run the indexer.
        # TODO: support both indexing modes
        RealTimeLogger.get().info("Indexing {}".format(graph_filename))
        subprocess.check_call(["{}/vg".format(bin_dir), "index", "-s", "-k",
            str(options.kmer_size), "-e", str(options.edge_max),
            "-t", str(job.cores), graph_filename])
            
        # Now save the indexed graph directory to the file store. It can be
        # cleaned up since only our children use it.
        index_dir_id = write_global_directory(job.fileStore, graph_dir,
            cleanup=True)
            
        # Add a child to actually save the graph to the output. Hack our own job
        # so that the actual alignment targets get added as a child of this, so
        # they happen after. TODO: massive hack!
        job.addChildJobFn(save_indexed_graph, options, index_dir_id,
            index_key, cores=1, memory="10G", disk="50G")
            
    RealTimeLogger.get().info("Done making children")
                    
    for sample in samples_to_run:
        # Split out over each sample that needs to be run
        
        # For each sample, know the FQ name
        sample_fastq = "{}/{}/{}.bam.fq".format(region_dir, sample, sample)
        
        # And know where we're going to put the output
        alignment_file_key = "{}/{}.gam".format(alignment_dir, sample)
        stats_file_key = "{}/{}.json".format(stats_dir, sample)
        
        RealTimeLogger.get().info("Queueing alignment of {} to {} {}".format(
            sample, graph_name, region))
    
        # Go and bang that input fastq against the correct indexed graph.
        # Its output will go to the right place in the output store.
        job.addChildJobFn(run_alignment, options, bin_dir_id, region,
            index_dir_id, sample_fastq, alignment_file_key, stats_file_key,
            cores=16, memory="100G", disk="50G")
            
def save_indexed_graph(job, options, index_dir_id, output_key):
    """
    Save the index dir tar file in the given output key.
    
    Runs as a child to ensure that the global file store can actually
    produce the file when asked (because within the same job, depending on Toil
    guarantees, it might still be uploading).
    
    """
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    sample_store = IOStore.get(options.sample_store)
    out_store = IOStore.get(options.out_store)
    
    # Get the tar.gz file
    local_path = job.fileStore.readGlobalFile(index_dir_id)
    
    # Save it as output
    out_store.write_output_file(local_path, output_key)
    
   
def run_alignment(job, options, bin_dir_id, region, index_dir_id,
    sample_fastq_key, alignment_file_key, stats_file_key):
    """
    Align the the given fastq from the input store against the given indexed
    graph (in the file store as a directory) and put the GAM and statistics in
    the given output keys in the output store.
    
    """
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    sample_store = IOStore.get(options.sample_store)
    out_store = IOStore.get(options.out_store)
    
    # Download the binaries
    bin_dir = "{}/bin".format(job.fileStore.getLocalTempDir())
    read_global_directory(job.fileStore, bin_dir_id, bin_dir)
    
    # Download the indexed graph to a directory we can use
    graph_dir = "{}/graph".format(job.fileStore.getLocalTempDir())
    read_global_directory(job.fileStore, index_dir_id, graph_dir)
    
    # We know what the vg file in there will be named
    graph_file = "{}/graph.vg".format(graph_dir)
    
    # Also we need the sample fastq
    fastq_file = "{}/input.fq".format(job.fileStore.getLocalTempDir())
    sample_store.read_input_file(sample_fastq_key, fastq_file)
    
    # And temp files for our aligner output and stats
    output_file = "{}/output.gam".format(job.fileStore.getLocalTempDir())
    stats_file = "{}/stats.json".format(job.fileStore.getLocalTempDir())
    
    # How long did the alignment take to run, in seconds?
    run_time = None
    
    # Open the file stream for writing
    with open(output_file, "w") as alignment_file:
    
        # Start the aligner and have it write to the file
        
        # Plan out what to run
        vg_parts = ["{}/vg".format(bin_dir), "map", "-f", fastq_file, "-i",
            "-n3", "-M2", "-t", str(job.cores), "-k", str(options.kmer_size),
            graph_file]
        
        RealTimeLogger.get().info("Running VG: {}".format(" ".join(vg_parts)))
        
        # Mark when we start the alignment
        start_time = timeit.default_timer()
        process = subprocess.Popen(vg_parts, stdout=alignment_file)
            
        if process.wait() != 0:
            # Complain if vg dies
            raise RuntimeError("vg died with error {}".format(
                process.returncode))
                
        # Mark when it's done
        end_time = timeit.default_timer()
        run_time = end_time - start_time
        
                
    RealTimeLogger.get().info("Aligned {}".format(output_file))
           
    # Read the alignments in in JSON-line format
    view = subprocess.Popen(["{}/vg".format(bin_dir), "view", "-aj",
        output_file], stdout=subprocess.PIPE)
       
    # Count up the stats: total reads, total mapped at all, total multimapped,
    # primary alignment score counts, secondary alignment score counts, and
    # aligner run time in seconds.
    
    stats = {
        "total_reads": 0,
        "total_mapped": 0,
        "total_multimapped": 0,
        "primary_scores": collections.Counter(),
        "primary_mismatches": collections.Counter(),
        "secondary_scores": collections.Counter(),
        "secondary_mismatches": collections.Counter(),
        "run_time": run_time
    }
        
    last_alignment = None
        
    for line in view.stdout:
        # Parse the alignment JSON
        alignment = json.loads(line)
        
        if alignment.has_key("score"):
            # This alignment is aligned.
            # Grab its score
            score = alignment["score"]
        
            # Calculate the mismatches
            length = len(alignment["sequence"])
            matches = 0
            for mapping in alignment.get("path", {}).get("mapping", []):
                for edit in mapping.get("edit", []):
                    if (not edit.has_key("sequence") and 
                        edit.get("to_length", None) == edit.get(
                        "from_length", None)):
                        
                        # We found a perfect match edit. Grab its length
                        matches += edit["from_length"]
                        
            # Calculate mismatches as what's left
            mismatches = length - matches
                    
        
            if alignment.get("is_secondary", False):
                # It's a multimapping. We can have max 1 per read, so it's a
                # multimapped read.
                
                if (last_alignment is None or 
                    last_alignment.get("name") != alignment.get("name") or 
                    last_alignment.get("is_secondary", False)):
                
                    # This is a secondary alignment without a corresponding primary
                    # alignment (which would have to be right before it given the
                    # way vg dumps buffers
                    raise RuntimeError("{} secondary alignment comes after "
                        "alignment of {} instead of corresponding primary "
                        "alignment\n".format(alignment.get("name"), 
                        last_alignment.get("name") if last_alignment is not None 
                        else "nothing"))
                
                # Log its stats as multimapped
                stats["total_multimapped"] += 1
                stats["secondary_scores"][score] += 1
                stats["secondary_mismatches"][mismatches] += 1
            else:
                # Log its stats as primary. We'll get exactly one of these per
                # read with any mappings.
                stats["total_mapped"] += 1
                stats["primary_scores"][score] += 1
                stats["primary_mismatches"][mismatches] += 1
                
                # We won't see an unaligned primary alignment for this read, so
                # count the read
                stats["total_reads"] += 1
        
        elif not alignment.get("is_secondary", False):
            # We have an unmapped primary "alignment"
            
            # Count the read by its primary alignment
            stats["total_reads"] += 1
            
        # Save the alignment for checking for wayward secondaries
        last_alignment = alignment
                
    with open(stats_file, "w") as stats_handle:
        # Save the stats as JSON
        json.dump(stats, stats_handle)
        
    # Now send the output files (alignment and stats) to the output store where
    # they belong.
    out_store.write_output_file(output_file, alignment_file_key)
    out_store.write_output_file(stats_file, stats_file_key)
    
        
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    if len(args) == 2 and args[1] == "--test":
        # Run the tests
        return doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    RealTimeLogger.start_master()
    
    # Pre-read the input file so we don't try to send file handles over the
    # network.
    options.server_list = list(options.server_list)
    
    # Make a root job
    root_job = Job.wrapJobFn(run_all_alignments, options,
        cores=1, memory="4G", disk="50G")
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print("All jobs completed successfully")
    
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

