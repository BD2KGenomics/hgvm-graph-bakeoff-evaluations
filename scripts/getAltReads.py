#!/usr/bin/env python2.7
"""
getAltReads.py: download reads for regions of interest from the 1000 Genomes FTP
server in parallel, using Toil.

"""

import argparse, sys, os, os.path, random, collections, shutil, itertools, glob
import urllib2, urlparse, ftplib, fnmatch, subprocess
import json, logging, logging.handlers, SocketServer, struct, socket, threading
import time, select

from toil.job import Job
import tsv

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
    parser.add_argument("--reference_metadata", 
        default="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
        "GCA_000001405.17_GRCh38.p2/"
        "GCA_000001405.17_GRCh38.p2_assembly_structure/"
        "all_alt_scaffold_placement.txt",
        help="URL to download the reference metadata from")
    parser.add_argument("--regions", nargs="*", 
        default=["BRCA1", "BRCA2", "CENX", "MHC", "SMA", "LRC_KIR"],
        help="region names to download reads for")
    parser.add_argument("--sample_ftp_root",
        default=("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
        "1000_genomes_project/data"), 
        help="FTP directory to scan for samples")
    parser.add_argument("--population_pattern", default="*", 
        help="fnmatch-style pattern for population names")
    parser.add_argument("--sample_pattern", default="NA*", 
        help="fnmatch-style pattern for sample names")
    parser.add_argument("--file_pattern", default="*.low_coverage.cram", 
        help="fnmatch-style pattern for read files in sample directories")
    parser.add_argument("--index_suffix", default=".crai", 
        help="suffix to add to sample files to get index files")
    parser.add_argument("--min_indexed_contigs", type=int, default=1000,
        help="reject samples with indexes not covering this many contigs")
    parser.add_argument("--sample_limit", type=int, default=float("inf"), 
        help="number of matching samples to download")
    parser.add_argument("--ftp_retry", type=int, default=float("inf"), 
        help="number of times to retry sample downloads")
    parser.add_argument("out_dir",
        help="output directory to create and fill with per-region BAM files")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)


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
    
    # We keep the master host and port as class data
    master_host = None
    master_port = None
    
    # Also the logger
    logger = None
  
    @classmethod
    def set_master(cls, options):
        """
        Set the master info.
        
        """
        
        cls.master_host = options.log_host
        cls.master_port = options.log_port
        
        
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
            cls.logger.addHandler(JSONDatagramHandler(cls.master_host,
                cls.master_port))
        
        return cls.logger

def backoff_times(retries=float("inf"), base_delay=300):
    """
    A generator that yields times for exponential back-off. Always yields 0
    first, and then if you have nonzero retries it yields exponentially but
    randomly increasing times in seconds to wait before trying again, stopping
    at the specified number of retries (and continuing forever by default).
    
    You have to do the error catching and sleeping yourself.
    
    Raises an exception if you use up all your backoff times.
    """
    
    # Don't wait at all before the first try
    yield 0
    
    # What retry are we on?
    try_number = 1
    
    # Make a delay that increases
    delay = float(base_delay) * 2
    
    while try_number <= retries:
        # Wait a random amount between 0 and 2^try_number * base_delay
        yield random.uniform(base_delay, delay)
        delay *= 2
        try_number += 1
        
    raise RuntimeError("Ran out of retries")

def ftp_connect(url, retries=float("inf")):
    """
    Connect to an FTP server and go to the specified directory with FTPlib.
    
    Return the ftplib connection and the path.
    
    On errors like timeouts, retry up to the given number of retries (infinite
    by default) with exponential back-off.
    
    """
    
    for delay in backoff_times(retries=retries):
        if delay > 0:
            # We have to wait before trying again
            RealTimeLogger.get().info("Retry after {} seconds".format(delay))
            time.sleep(delay)
        try:
        
            ftp_info = urlparse.urlparse(url)
            assert(ftp_info.scheme == "ftp")
            
            # Connect to the server
            ftp = ftplib.FTP(ftp_info.netloc, ftp_info.username, ftp_info.password)
            
            # Log in
            ftp.login()
                
            # Go to the right directory
            ftp.cwd(ftp_info.path)
            
            return ftp, ftp_info.path
            
        except IOError as e:
            # Something went wrong doing the IO
            RealTimeLogger.get().warning(
                "Retry after FTP setup IO error: {}".format(e))
            

def explore_path(ftp, path, pattern):
    """
    Using the given FTP server connection, explore the given path recursively,
    looking for all files that match the given fnmatch pattern. Yields the paths
    of those files.
    
    """
    
    try:
        ftp.cwd(path)
        
        # List of everything to recurse into
        to_recurse_on = []
        
        try:
            listing = ftp.nlst()
        except ftplib.error_proto:
            # For some reason this built-in function sometimes doesn't work. Try
            # doing it ourselves.
            
            RealTimeLogger.get().warning("FTPlib nlst on {} failed. Clearing "
                "and retrying.".format(path))
            
            while True:
                # Clear out any crap that may be coming over the ftp connection.
                # Use a 2 second timeout
                flags = select.select([ftp.sock], [], [], 2)
                
                if not flags[0]:
                    # Nothing came
                    break
                    
                # Otherwise read some data and select again
                got = ftp.sock.recv(1024)
                RealTimeLogger.get().warning("Extra data: {}".format(got))
            
            
            # Now manually do the listing
            listing = []
            ftp.retrlines("NLST", lambda line: listing.append(line))
        
        for subitem in listing:
            # We don't know if these are files or directories.
            
            subitem_path = "{}/{}".format(path, subitem)
            
            print(subitem_path, pattern)
            
            if fnmatch.fnmatchcase(subitem, pattern):
                # This is a matching thing!
                print("Match!")
                yield subitem_path
                
            # Recurse on everything, even things that match the pattern, in case
            # they are directories.
            to_recurse_on.append(subitem_path)
            
        for item in to_recurse_on:
            # Recurse down on each thing in turn
            for found in explore_path(ftp, item, pattern):
                # Yield anything we find there
                yield found
            # Then return to this directory so we're in a known place
            ftp.cwd(path)
            
    except ftplib.error_perm as e:
        error_code = int(e.args[0][:3])
        if error_code == 550:
            # We expect to do a lot of CWD-ing to files, raising this
            pass
        else:
            raise e
        
def count_indexed_contigs(index_url, retries):
    """
    Given the URL of a .crai index file, count the number of distinct contigs in
    the index and return it.
    
    Basically: curl url | zcat | cut -f1 | uniq | wc -l
    """
    
    for delay in backoff_times(retries=retries):
        if delay > 0:
            # We have to wait before trying again
            RealTimeLogger.get().info("Retry after {} seconds".format(delay))
            time.sleep(delay)
        try:
        
            # Hold all the popen objects we need for this
            tasks = []
            
            # Do the download
            tasks.append(subprocess.Popen(["curl", "--silent", index_url],
                stdout=subprocess.PIPE))
            
            # Pipe through zcat
            tasks.append(subprocess.Popen(["zcat"], stdin=tasks[-1].stdout,
                stdout=subprocess.PIPE))
            
            # And cut
            tasks.append(subprocess.Popen(["cut", "-f1"], stdin=tasks[-1].stdout,
                stdout=subprocess.PIPE))
                
            # And uniq
            tasks.append(subprocess.Popen(["uniq"], stdin=tasks[-1].stdout,
                stdout=subprocess.PIPE))
                
            # And wc
            output = subprocess.check_output(["wc", "-l"],
                stdin=tasks[-1].stdout)
            
            # Did we make it through all the tasks OK?
            success = True
            for task in tasks:
                if task.wait() != 0:
                    RealTimeLogger.get().warning(
                        "Pipeline failed: {}".format(index_url)) 
                    success = False
                    break
                    
            if success:
                # The pipeline worked.
                # Parse and return the line count
                return int(output.strip())
            
        except subprocess.CalledProcessError as e:
            # Something went wrong doing the IO
            RealTimeLogger.get().warning(
                "Pipeline failed: {}".format(e)) 
    
def downloadAllReads(job, options):
    """
    Download all the reads for the regions.
    
    """
    
    # Initialize logging
    RealTimeLogger.set_master(options)
    
    RealTimeLogger.get().info("Starting download")
    
    # First make the output directory
    if not os.path.exists(options.out_dir):
        try:
            # Make it if it doesn't exist
            os.makedirs(options.out_dir)
        except OSError:
            # If you can't make it, maybe someone else did?
            pass
            
    # Whatever happens, it needs to exist here
    assert(os.path.exists(options.out_dir) and os.path.isdir(options.out_dir))
    
    # Holds the chromosome number for each region?
    region_chromosomes = {}
    # Holds the minimum start position for each region on its chromosome
    region_starts = collections.defaultdict(lambda: float("inf"))
    # Holds the maximum stop position for each region on its chromosome
    region_stops = collections.defaultdict(lambda: float("-inf"))
    
    # Holds the contig:start-end string for each alt in a given region
    # The reference range gets added in last
    ranges_by_region = collections.defaultdict(list)
    
    # Hard-code some regions that aren't real
    ranges_by_region["BRCA1"] = ["chr17:43044294-43125482"]
    ranges_by_region["BRCA2"] = ["chr13:32314861-32399849"]
    ranges_by_region["CENX"] = ["chrX:58605580-62412542"]
    
    # Read the reference database
    database = tsv.TsvReader(urllib2.urlopen(options.reference_metadata))
    
    for parts in database:
        # Parse out all the info for this alt and its parent chromosome
        region_name = parts[7]
        # Grab the chromosome ("1" or "X") that's the parent
        parent_chromosome = parts[5]
        parent_start = int(parts[11])
        parent_stop = int(parts[12])
        alt_contig = parts[3]
        alt_start = int(parts[9])
        alt_stop = int(parts[10])
        
        # Note the region start, stop, and parent chromosome number
        region_chromosomes[region_name] = parent_chromosome
        region_starts[region_name] = min(region_starts[region_name],
            parent_start)
        region_stops[region_name] = max(region_stops[region_name],
            parent_stop)
            
        # Turn the alt name into the proper format (GL000251.2 to
        # chr6_GL000251v2_alt)
        name_parts=alt_contig.split(".")
        fixed_alt_contig = "chr{}_{}v{}_alt".format(parent_chromosome,
            name_parts[0], name_parts[1])
            
        # Add it to the list for its region
        ranges_by_region[region_name].append("{}:{}-{}".format(
            fixed_alt_contig, alt_start, alt_stop))
                
    for region_name in region_chromosomes.iterkeys():
        # Add in the reference ranges that all the alts are alternatives to
        ranges_by_region[region_name].append("chr{}:{}-{}".format(
            region_chromosomes[region_name], region_starts[region_name], 
            region_stops[region_name]))
            

    ftp, root_path = ftp_connect(options.sample_ftp_root)
    
    # Calculate the FTP base URL (without directory). We need it later for
    # turning found index files into URLs.
    base_url = options.sample_ftp_root[:-len(root_path)]
    
    # Dump the good data files for samples
    good_samples = open("{}/good.txt".format(options.out_dir), "w")
    
    # Grab all the population names that match the population pattern
    population_names = [n for n in ftp.nlst() if fnmatch.fnmatchcase(n,
        options.population_pattern)]
        
    # TODO: We'll go through them in this order, so if you want a representative
    # subsampling, add some shuffle here or something.

    for population_name in population_names:
    
        # For each of those, we need to get samples
        # Go to the population root
        ftp.cwd("{}/{}".format(root_path, population_name))
        
        # Grab all the sample names that match the sample name pattern.
        # Hopefully there aren't too many.
        sample_names = [n for n in ftp.nlst() if fnmatch.fnmatchcase(n,
            options.sample_pattern)]
        
        # This holds URLs to data files (BAM/CRAM) with indexes that are on a
        # sufficient number of contigs, by sample name. We take the first
        # sufficiently good file for any sample.
        sample_file_urls = {}
        
        # TODO: handle failures during explore_path?
        for sample_name in sample_names:
            # For every sample
            
            print("Try {}".format(sample_name))
            
            for data_name in explore_path(ftp, "{}/{}/{}".format(root_path,
                population_name, sample_name), options.file_pattern):
                # Find its data files (there may be several)
                
                # Get the index for each
                index_name = data_name + options.index_suffix
                
                print(index_name)
                
                # Count up the contigs it indexes over
                indexed_contigs = count_indexed_contigs("{}/{}".format(
                    base_url, index_name), options.ftp_retry)
                    
                if indexed_contigs >= options.min_indexed_contigs:
                    # This file for this sample is good enough
                    sample_file_urls[sample_name] = "{}/{}".format(base_url,
                        data_name)
                    
                    RealTimeLogger.get().info(
                        "Sample {} has index of {} contigs".format(sample_name,
                        indexed_contigs))
                    # Add the sample to the file we spit out
                    good_samples.write("{}\n".format(sample_name))
                    
                    # Don't finish exploring the path
                    break
                    
                else:
                    # Complain
                    RealTimeLogger.get().warning(
                        "Sample {} has index on too few contigs ({})."
                        "Skipping!".format(sample_name, indexed_contigs))
                    
            if len(sample_file_urls) >= options.sample_limit:
                # We got enough. Don't finish this population
                break
                
        if len(sample_file_urls) >= options.sample_limit:
                # We got enough. Don't check more populations
                break
            
    good_samples.close()
            
    RealTimeLogger.get().info("Got {} sample URLs".format(
        len(sample_file_urls)))
    
    if(options.sample_limit < float("inf")):
        # Make sure we got as many as we wanted.
        assert(len(sample_file_urls) == options.sample_limit)
    
    for region_name in options.regions:
        for sample_name, sample_url in sample_file_urls.iteritems():
        
            # Make sure the sample directory exists
            sample_dir = "{}/{}/{}".format(options.out_dir, region_name,
                sample_name)
            
            if not os.path.exists(sample_dir):
                try:
                    # Make it if it doesn't exist
                    os.makedirs(sample_dir)
                except OSError:
                    # If you can't make it, maybe someone else did?
                    pass
                    
            assert(os.path.exists(sample_dir) and os.path.isdir(sample_dir))
            
            # Where will this sample's BAM for this region go?
            bam_filename = "{}/{}.bam".format(sample_dir, sample_name)
            
            RealTimeLogger.get().info("Making child for {} x {}: {}".format(
                region_name, sample_name, sample_url))
            
            # Now kick off a job to download all the ranges for the region in
            # parallel for this sample, and then concatenate them together. Tell
            # it to save the results to a file on a shared filesystem.
            job.addChildJobFn(downloadRegion, options, region_name, 
                sample_url, ranges_by_region[region_name], bam_filename, 
                cores=1, memory="1G", disk=0)
                
    RealTimeLogger.get().info("Done making children")
   
def downloadRegion(job, options, region_name, file_url, range_list,
    bam_filename):
    """
    Download all the ranges given for the given region from the given sample
    data file URL, concatenate them, and save them to the given BAM.
    
    """
    
    RealTimeLogger.set_master(options)
    
    # TODO: implement
    RealTimeLogger.get().info("Supposed to download {} ranges {} to {}".format(
        file_url, range_list, bam_filename))
        
    # For each region, we get a promise of a returned file store file ID. This
    # holds them.
    part_promises = []
    
    # We just have the one file
    
    for range_string in range_list:
        # For every range we want from it...
        
        # Set up a child job to grab it that returns a file store ID for the
        # BAM file it gets. Right now these are promises, but they get
        # filled in later.
        part_promises.append(job.addChildJobFn(downloadRange, options,
            file_url, range_string, cores=1, memory="1G",
            disk="50G").rv())
                
    # Make a follow-on that concatenates the parts together
    job.addFollowOnJobFn(concatAndSortBams, options, part_promises,
        bam_filename, cores=1, memory="4G", disk="50G")
        
        
def downloadRange(job, options, file_url, range_string):
    """
    Download the given range from the given HTSlib file, put the resultiung BAM
    in the file store, and return its file ID.
    
    """
    
    RealTimeLogger.set_master(options)
    
    # Where should we save the bam locally?
    bam_filename = "{}/download.bam".format(job.fileStore.getLocalTempDir())
    
    for delay in backoff_times(retries=options.ftp_retry):
        if delay > 0:
            # We have to wait before trying again
            RealTimeLogger.get().info("Retry after {} seconds".format(delay))
            time.sleep(delay)
        try:
            # Try running the download
            RealTimeLogger.get().info("Trying to download {} from {}".format(
                range_string, file_url))
            subprocess.check_call(["samtools", "view", "-b", "-o", bam_filename,
                file_url, range_string])
                
            # Make sure there's actually reads in its file.
            lines = subprocess.check_output(["samtools", "flagstat",
                bam_filename])
                
            # Badly parse the flagstat output for the very first number (QC-
            # passed reads)
            reads_found = int(lines.split()[0])
                
            if reads_found > 0:
                
                # If we get here it worked
                break
            else:
                # Complain we downloaded and there were no reads
                RealTimeLogger.get().warning(
                    "No reads! Need to retry download of {}".format(file_url))
                    
            # TODO: Catch a partial download somehow
        except subprocess.CalledProcessError as e:
                # Complain we need to retry
                RealTimeLogger.get().warning(
                    "Need to retry download of {}".format(file_url))
                # But keep looping
                    
    # Put the BAM in the file store with a new ID
    file_id = job.fileStore.writeGlobalFile(bam_filename, cleanup=False)

    RealTimeLogger.get().info("Wrote file to ID {}".format(file_id))
    
    return file_id
    
def concatAndSortBams(job, options, bam_ids, output_filename):
    """
    Takes in a list of BAM file IDs in the file store, concatenates and sorts
    them by template name, and puts the result in the specified shared
    filesystem file.
    
    Also creates an interleaved FASTQ file of the paired reads at
    <output_filename>.fq
    
    """
    
    RealTimeLogger.set_master(options)
    
    RealTimeLogger.get().info("Building {} from IDs {}".format(output_filename,
        bam_ids))
    
    # Where do we put the concatenated bam?
    concat_filename = "{}/concat.bam".format(
        job.fileStore.getLocalTempDir())
    
    if len(bam_ids) > 1:
        RealTimeLogger.get().info("Concatenating...")
        
        # What input BAMs do we want?
        input_files = [job.fileStore.readGlobalFile(bam_id)
            for bam_id in bam_ids]
        
        # Do the concatenation
        subprocess.check_call(["samtools", "cat", "-o", concat_filename] +
            input_files)
    elif len(bam_ids) == 1:
        RealTimeLogger.get().info("Skip concatenation...")
        
        # Just grab the file for the one BAM
        job.fileStore.readGlobalFile(bam_ids[0], concat_filename)
    else:
        raise RuntimeError("Got no BAM parts!")
    
    RealTimeLogger.get().info("Creating {}".format(output_filename))
    
    # Make up a prefix for sort. It gets .bam appended
    sort_prefix = "{}/sort".format(job.fileStore.getLocalTempDir())
    
    # Do the sort into the temp file
    subprocess.check_call(["samtools", "sort", "-n", concat_filename,
        sort_prefix])
    
    
    
    # Convert to FASTQ
    # Decide on the temp filename
    fastq_filename = "{}/reads.fq".format(job.fileStore.getLocalTempDir())
    RealTimeLogger.get().info("Creating {}.fq".format(output_filename))
    
    # View the sam and pipe it through the sort of smart converter
    view = subprocess.Popen(["samtools", "view", "{}.bam".format(sort_prefix)],
        stdout=subprocess.PIPE)
    subprocess.check_call(["./smartSam2Fastq.py", "--interleaved", "--fq1",
        fastq_filename], stdin=view.stdout)
    
    # Wait and detect errors
    view.wait()
    assert(view.returncode == 0)
    
    # Move it
    shutil.move("{}.bam".format(sort_prefix), output_filename)
    shutil.move(fastq_filename, "{}.fq".format(output_filename))
        
def main():
    options = parse_args(sys.argv) # This holds the nicely-parsed options object
    
    logging.basicConfig(level=logging.DEBUG)
    
    # Start up the logging server
    logging_server = SocketServer.ThreadingUDPServer(("0.0.0.0", 0),
        LoggingDatagramHandler)
        
    print("Starting server")
    # Set up a thread to do all the serving in the background and exit when we
    # do
    server_thread = threading.Thread(target=logging_server.serve_forever)
    server_thread.daemon = True
    server_thread.start()
    print("Server started")
    
    # HACK: Set options for logging and pass the options to every target
    options.log_host = socket.getfqdn()
    options.log_port = logging_server.server_address[1]
    
    RealTimeLogger.set_master(options)
    
    logger = RealTimeLogger.get()
    
    # Make the root job
    root_job = Job.wrapJobFn(downloadAllReads, options, 
        cores=1, memory="1G", disk=0)
        
    print("Sending log from master")
    logger.info("This is the master")
    
    # Run Toil
    Job.Runner.startToil(root_job,  options)
        
    logging_server.shutdown()
    server_thread.join()

if __name__=="__main__":
    sys.exit(main())
