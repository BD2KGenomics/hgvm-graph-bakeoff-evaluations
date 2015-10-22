#!/usr/bin/env python2.7
"""
histogram: plot a histogram of a file of numbers. Numbers can be floats, one per
line. Lines with two numbers are interpreted as pre-counted, with the number of
repeats of the first being given by the second.

Re-uses sample code and documentation from 
<http://users.soe.ucsc.edu/~karplus/bme205/f12/Scaffold.html>
"""

import argparse, sys, os, itertools, math, numpy, collections
import matplotlib, matplotlib.ticker

def intify(x):
    """
    Turn an integral float into an int, if applicable.
    """
    
    if isinstance(x, float) and x.is_integer():
        return int(x)
    return x
    
def draw_labels(bin_counts, bar_patches, size=None):
    """
    Put the given count labels on the given bar patches, on the current axes.
    Takes an optional font size.
    
    """
    
    from matplotlib import pyplot
    
    # Grab the axes
    axes = pyplot.gca()

    for bin_count, bar_patch in itertools.izip(bin_counts, bar_patches):
        
        if(bin_count.is_integer()):
            # Intify if applicable
            bin_count = int(bin_count)
        
        # Label each bar
        if bin_count == 0:
            # Except those for empty bins
            continue
            
        # Find the center of the bar
        bar_center_x = bar_patch.get_x() + bar_patch.get_width() / float(2)
        # And its height
        bar_height = bar_patch.get_height()
        
        # Label the bar
        axes.annotate("{:,}".format(bin_count), (bar_center_x, bar_height), 
            ha="center", va="bottom", rotation=45, xytext=(0, 5), 
            textcoords="offset points", size=size)

def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Now add all the options to it
    parser.add_argument("data", type=argparse.FileType('r'), nargs="+",
        help="the file to read")
    parser.add_argument("--redPortion", type=float, action="append", default=[],
        help="portion of each bin to color red")
    parser.add_argument("--redWeight", type=float, action="append", default=[],
        help="value to plot in red in each bin")
    parser.add_argument("--title", default="Histogram",
        help="the plot title")
    parser.add_argument("--x_label", default="Value",
        help="the plot title")
    parser.add_argument("--y_label", default="Number of Items (count)",
        help="the plot title")
    parser.add_argument("--bins", type=int, default=10,
        help="the number of histogram bins")
    parser.add_argument("--x_min", "--min", type=float, default=None,
        help="minimum value allowed")
    parser.add_argument("--x_max", "--max", type=float, default=None,
        help="maximum value allowed")
    parser.add_argument("--y_min", type=int, default=None,
        help="minimum count on plot")
    parser.add_argument("--y_max", type=int, default=None,
        help="maximum count on plot")
    parser.add_argument("--cutoff", type=float, default=None,
        help="note portion above and below a value")
    parser.add_argument("--font_size", type=int, default=12,
        help="the font size for text")
    parser.add_argument("--cumulative", action="store_true",
        help="plot cumulatively")
    parser.add_argument("--log", action="store_true",
        help="take the base-10 logarithm of values before plotting histogram")
    parser.add_argument("--log_counts", "--logCounts", action="store_true",
        help="take the logarithm of counts before plotting histogram")
    parser.add_argument("--stats", action="store_true",
        help="print data stats")
    parser.add_argument("--save",
        help="save figure to the given filename instead of showing it")
    parser.add_argument("--dpi", type=int, default=300,
        help="save the figure with the specified DPI, if applicable")
    parser.add_argument("--sparse_ticks", action="store_true",
        help="use sparse tick marks")
    parser.add_argument("--label", action="store_true",
        help="label bins with counts")
    parser.add_argument("--label_size", type=float,
        help="bin count label font size")
    parser.add_argument("--no_n", dest="show_n", action="store_false",
        help="don't add n value to title")
    parser.add_argument("--normalize", action="store_true",
        help="normalize to total weight of 1")
    parser.add_argument("--line", action="store_true",
        help="draw a line instead of a barchart")
    parser.add_argument("--points", action="store_true",
        help="draw points instead of a barchart")
    parser.add_argument("--labels", nargs="*", default=[],
        help="labels for each data file histogram")
    parser.add_argument("--width", type=float, default=8,
        help="plot width in inches")
    parser.add_argument("--height", type=float, default=6,
        help="plot height in inches")
    
        
    return parser.parse_args(args)
    
def filter2(criterion, key_list, other_list):
    """
    Filter two lists of corresponding items based on some function of the first
    list.
    
    """
    
    # Make the output lists
    out1 = []
    out2 = []
    
    for key_val, other_val in itertools.izip(key_list, other_list):
        # Pair up the items
        if criterion(key_val):
            # The key passed the filter, so take both.
            out1.append(key_val)
            out2.append(other_val)
            
    return out1, out2

def main(args):
    """
    Parses command line arguments, and plots a histogram.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    if options.save is not None:
        # Set up plot for use in headless mode if we just want to save. See
        # <http://stackoverflow.com/a/2766194/402891>. We need to do this before
        # we grab pyplot.
        matplotlib.use('Agg')
        
    from matplotlib import pyplot
    
    # Make the figure with the appropriate size and DPI.
    pyplot.figure(figsize=(options.width, options.height), dpi=options.dpi)
    
    # This will hold a list of lists of data, weight pairs
    all_data = []
    
    for data_file in options.data:
        # This holds all the data values
        data = []
        
        # This holds a weight for each to allow for pre-counted values
        weights = []
        
        for line_number, line in enumerate(data_file):
            # Split each line
            parts = line.split()
            
            if len(parts) == 1:
                # This is one instance of a value
                data.append(float(parts[0]))
                weights.append(1)
            elif len(parts) == 2:
                # This is multiple instances of a value
                data.append(float(parts[0]))
                weights.append(float(parts[1]))
            else:
                raise Exception("Wrong number of fields on line {}".format(
                    line_number + 1))
        
        # TODO: Make this streaming instead of loading.
        
        # Strip NaNs and Infs
        data, weights = filter2(lambda x: x < float("+inf") and 
            x > float("-inf"), data, weights)
            
        # Strip out data that has weight 0
        weights, data = filter2(lambda x: x > 0, weights, data)
        
        # Save the lists
        all_data.append((data, weights))
        
    
    # Calculate our own bins, over all the data. First we need the largest and
    # smallest observed values. The fors in the comprehension have to be in
    # normal for loop order and not the other order.
    bin_min = options.x_min if options.x_min is not None else min((item
        for pair in all_data for item in pair[0]))
    bin_max = options.x_max if options.x_max is not None else max((item
        for pair in all_data for item in pair[0]))
    
    if options.log:
        # Do our bins in log space, so they look evenly spaced on the plot.
        bin_max = math.log10(bin_max)
        bin_min = math.log10(bin_min)
    
    # Work out what step we should use between bin edges
    bin_step = (bin_max - bin_min) / float(options.bins)
    # Work out where the bin edges should be
    bins = [bin_min + bin_step * i for i in xrange(options.bins + 1)]
    # Work out where the bin centers should be
    bin_centers = [left_edge + bin_step / 2.0 for left_edge in bins[:-1]]    
    
    if options.log:
        # Bring bins back into data space
        bins = [math.pow(10, x) for x in bins]
        bin_centers = [math.pow(10, x) for x in bin_centers]
    
    for (data, weights), label, line_style in itertools.izip(all_data,
        itertools.chain(options.labels, itertools.repeat(None)), 
        itertools.cycle(['-', '--', ':', '-.'])):
        
        # For each set of data and weights that we want to plot, and the label
        # it needs (or None)...
            
        # Apply the limits
        if options.x_min is not None:
            data, weights = filter2(lambda x: x > options.x_min, data, weights)
        if options.x_max is not None:
            data, weights = filter2(lambda x: x < options.x_max, data, weights)
            
        # Let's condense down by summing all weights for values
        total_weight = collections.defaultdict(lambda: 0)
        # We need a float here so we don't get int division later.
        total_weight_overall = float(0)
        
        for value, weight in itertools.izip(data, weights):
            # Sum up the weights for each value
            total_weight[value] += weight
            # And overall
            total_weight_overall += weight
        
        # Unpack the data and summed weights
        data, weights = zip(*total_weight.items())
        
        if options.normalize and total_weight_overall > 0:
            # Normalize all the weight to 1.0 total weight.
            weights = [w / total_weight_overall for w in weights]
           
        # Work out how many samples there are
        samples = intify(sum(weights))
            
        if options.stats:
            # Compute and report some stats
            data_min = numpy.min(data)
            data_min_count = weights[numpy.argmin(data)]
            data_max = numpy.max(data)
            data_max_count = weights[numpy.argmax(data)]
            # The mode is the data item with maximal count
            data_mode = data[numpy.argmax(weights)]
            data_mode_count = numpy.max(weights)
            
            # Intify floats pretending to be ints
            data_min = intify(data_min)
            data_min_count = intify(data_min_count)
            data_max = intify(data_max)
            data_max_count = intify(data_max_count)
            data_mode = intify(data_mode)
            data_mode_count = intify(data_mode_count)
        
            # TODO: median, mean
            
            print("Min: {} occurs {} times".format(data_min, data_min_count))
            print("Mode: {} occurs {} times".format(data_mode, data_mode_count))
            print("Max: {} occurs {} times".format(data_max, data_max_count))
            
            if options.cutoff is not None:
                # Work out how much weight is above and below the cutoff
                above = 0
                below = 0
                
                for value, weight in itertools.izip(data, weights):
                    if value > options.cutoff:
                        above += weight
                    else:
                        below += weight
                
                # Report the results wrt the cutoff.
                print "{} above {}, {} below".format(
                    above / total_weight_overall, options.cutoff, 
                    below / total_weight_overall)
            
        if options.line or options.points:
            # Do histogram binning manually
            
            # Do the binning
            bin_values, _ = numpy.histogram(data, bins=bins, weights=weights)
            
            if options.cumulative:
                # Calculate cumulative weights for each data point
                bin_values = numpy.cumsum(bin_values)
            
            if options.line:
                # Do the plot as a line. Make sure we start and end at 0.    
                pyplot.plot([bins[0]] + list(bin_centers) + [bins[-1]],
                    [0] + list(bin_values) + [0], label=label,
                    linestyle=line_style)
            if options.points:
                # Do the plot as points.   
                pyplot.scatter(bin_centers, bin_values, label=label)
            
            if options.log_counts:
                # Log the Y axis
                pyplot.yscale('log')
        
        else:
            # Do the plot. Do cumulative, or logarithmic Y axis, optionally.
            # Keep the bin total counts and the bar patches.
            bin_counts, _, bar_patches = pyplot.hist(data, bins,
                cumulative=options.cumulative, log=options.log_counts,
                weights=weights, alpha=0.5 if len(options.data) > 1 else 1.0,
                label=label)
            
        
            
    if len(options.labels) > 0:
        pyplot.legend()
        
    if len(options.redPortion) > 0:
        # Plot a red histogram over that one, modified by redPortion.
        
        red_data = []
        red_weights = []
        
        for item, weight in itertools.izip(data, weights):
            # For each item, what bin is it in?
            bin_number = int(item / bin_step)
            
            if bin_number < len(options.redPortion):
                # We have a potentially nonzero scaling factor. Apply that.
                weight *= options.redPortion[bin_number]
                
                # Keep this item.
                red_data.append(item)
                red_weights.append(weight)
        
        # Plot the re-weighted data with the same bins, in red
        red_counts, _, red_patches = pyplot.hist(red_data, bins,
            cumulative=options.cumulative, log=options.log_counts,
            weights=red_weights, color='#FF9696', hatch='/'*6)
            
        if options.label:
            # Label all the red portion-based bars
            draw_labels(red_counts, red_patches, size=options.label_size)
            
        
            
    if len(options.redWeight) > 0:
        # Plot a red histogram over that one, modified by redPortion.
        
        # Grab an item in each bin
        items = bins[0:len(options.redWeight)]
        
        # Plot the re-weighted data with the same bins, in red
        red_counts, _, red_patches = pyplot.hist(items, bins,
            cumulative=options.cumulative, log=options.log_counts,
            weights=options.redWeight, color='#FF9696', hatch='/'*6)
            
        if options.label:
            # Label all the red weight-based bars
            draw_labels(red_counts, red_patches, size=options.label_size)
            
        
        
    # StackOverflow provides us with font sizing. See
    # <http://stackoverflow.com/q/3899980/402891>
    matplotlib.rcParams.update({"font.size": options.font_size})
    if options.show_n:
        # Add an n value to the title
        options.title += " (n = {:,})".format(samples)
    pyplot.title(options.title)
    pyplot.xlabel(options.x_label)
    pyplot.ylabel(options.y_label)
    
    if options.log:
        # Set the X axis to log mode
        pyplot.xscale('log')
    
    if options.x_min is not None:
        # Set only the lower x limit
        pyplot.xlim((options.x_min, pyplot.xlim()[1]))
    if options.x_max is not None:
        # Set only the upper x limit
        pyplot.xlim((pyplot.xlim()[0], options.x_max))
        
    if options.y_min is not None:
        # Set only the lower y limit
        pyplot.ylim((options.y_min, pyplot.ylim()[1]))
    elif options.log_counts:
        # Make sure the default lower Y limit is 1 on log plots.
        pyplot.ylim((1, pyplot.ylim()[1]))
    if options.y_max is not None:
        # Set only the upper y limit
        pyplot.ylim((pyplot.ylim()[0], options.y_max))
        
    if options.sparse_ticks:
        # Set up tickmarks to have only 2 per axis, at the ends
        pyplot.gca().xaxis.set_major_locator(
            matplotlib.ticker.FixedLocator(pyplot.xlim()))
        pyplot.gca().yaxis.set_major_locator(
            matplotlib.ticker.FixedLocator(pyplot.ylim()))
            
    # Make sure tick labels don't overlap. See
    # <http://stackoverflow.com/a/20599129/402891>
    pyplot.gca().tick_params(axis="x", pad=0.5 * options.font_size)

    if options.label:
        # Label all the normal bars
        draw_labels(bin_counts, bar_patches, size=options.label_size)
    
    # Make everything fit
    pyplot.tight_layout()
    
    if options.save is not None:
        # Save the figure to a file
        pyplot.savefig(options.save, dpi=options.dpi)
    else:
        # Show the figure to the user
        pyplot.show()
        
    return 0

if __name__ == "__main__" :
    sys.exit(main(sys.argv))

