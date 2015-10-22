#!/usr/bin/env python2.7
"""
scatter: plot a scatterplot of a file of numbers. Numbers should be floats, two 
per line.

Re-uses sample code and documentation from 
<http://users.soe.ucsc.edu/~karplus/bme205/f12/Scaffold.html>
"""

import argparse, sys, os, itertools, math, collections, random, re
import matplotlib, matplotlib.ticker, matplotlib.cm, numpy

# Implementation of "natural" sorting from
# <http://stackoverflow.com/a/5967539/402891>
def atoi(text):
    """
    Turn an int string into a number, but leave a non-int string alone.
    """
    return int(text) if text.isdigit() else text

def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    """
    return [atoi(c) for c in re.split('(\d+)', text)]

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
    parser.add_argument("data", type=argparse.FileType('r'),
        help="the file to read")
    parser.add_argument("--dotplot", action="store_true",
        help="use categories for the x axis")
    parser.add_argument("--title", default="Scatterplot",
        help="the plot title")
    parser.add_argument("--x_label", default="X",
        help="the X axis label")
    parser.add_argument("--log_x", action="store_true",
        help="log X axis")
    parser.add_argument("--y_label", default="Y",
        help="the Y axis label")
    parser.add_argument("--log_y", action="store_true",
        help="log Y axis")
    parser.add_argument("--font_size", type=int, default=12,
        help="the font size for text")
    parser.add_argument("--min_x", type=float, default=None,
        help="lower limit of X axis")
    parser.add_argument("--max_x", type=float, default=None,
        help="upper limit of X axis")
    parser.add_argument("--min_y", type=float, default=None,
        help="lower limit of Y axis")
    parser.add_argument("--max_y", type=float, default=None,
        help="upper limit of Y axis")
    parser.add_argument("--save",
        help="save figure to the given filename instead of showing it")
    parser.add_argument("--dpi", type=int, default=300,
        help="save the figure with the specified DPI, if applicable")
    parser.add_argument("--sparse_ticks", action="store_true",
        help="use sparse tick marks")
    parser.add_argument("--lines", action="store_true",
        help="connect points together")
    parser.add_argument("--tsv", action="store_true",
        help="use only tabs as separators in input file")
    parser.add_argument("--no_sort", dest="sort", action="store_false",
        help="do not sort categories from input file")
    parser.add_argument("--categories", nargs="+", default=None,
        help="categories to plot, in order")
    parser.add_argument("--legend_overlay", default=None,
        help="display the legend overlayed on the graph at this location")
    parser.add_argument("--colors", nargs="+", default=None,
        help="use the specified Matplotlib colors")
    parser.add_argument("--markers", nargs="+", default=None,
        help="use the specified Matplotlib markers")
    parser.add_argument("--no_n", dest="show_n", action="store_false",
        help="don't add n value to title")
    parser.add_argument("--width", type=float, default=8,
        help="plot width in inches")
    parser.add_argument("--height", type=float, default=6,
        help="plot height in inches")
    
        
    return parser.parse_args(args)
    

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
    
    # Make the figure with the appropriate size.
    pyplot.figure(figsize=(options.width, options.height))
    
    # This holds pairs of x, y lists for each data series.
    series = collections.defaultdict(lambda: [list(), list()])
    
    # This holds the order in which series were first encountered
    initial_series_order = collections.OrderedDict()
    
    # Should we use series or not?
    use_series = False
    
    if options.dotplot:
        # Dotplots always need series
        use_series = True
    
    for line in options.data:
        # Unpack the line, splitting on tabs only if requested
        parts = line.split("\t" if options.tsv else None)
        
        if options.dotplot:
            # We parse a two-column name/sample format
            series_name = parts[0]
            y_value = float(parts[1])
            # We fill in the x values later according to the order of the
            # series.
            x_value = None
        elif len(parts) == 3:
            # We have a series name. Pull that out and use it.
            series_name = parts[0]
            x_value = float(parts[1])
            y_value = float(parts[2])
            # We should be using series.
            use_series = True
        else:
            # Use the default series
            series_name = ""
            x_value = float(parts[0])
            y_value = float(parts[1])
            
        
        # Put each coordinate component in the appropriate list.
        series[series_name][0].append(x_value)
        series[series_name][1].append(y_value)
        
        # Note this series in the ordering if we haven't already
        initial_series_order[series_name] = None
        
    
    if use_series:
        
        if options.categories is not None:
            # Fix up the options so we don't ask for categories with no points.
            for i in xrange(len(options.categories)):
                if series.has_key(options.categories[i]):
                    # Leave this one
                    continue
                    
                # Otherwise none it out from everything
                options.categories[i] = None
                
                if options.colors is not None:
                    options.colors[i] = None
                
                if options.markers is not None:
                    options.markers[i] = None
             
            # Now filter out everrything we noned out       
            options.categories = [x for x in options.categories
                if x is not None]
                
            if options.colors is not None:
                options.colors = [x for x in options.colors
                if x is not None]
                
            if options.markers is not None:
                options.markers = [x for x in options.markers
                if x is not None]
        
        if options.colors is not None:
            # Use user-specified colors
            colors = options.colors
        else:
            # Make up colors for each series. We have 7.
            colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        
        # Cycle them in case there are not enough    
        series_colors = itertools.cycle(colors)

        
        if options.markers is not None:
            # Use user-specified markers
            markers = options.markers
        else:
            # Make up symbols for the marker list. We have 11, for 7*11
            # combinations before we repeat.
            markers = ['o', 'v', '^', '<', '>', 's', '+', 'x', 'D', '|', '_']
            # Make sure they're in a good order
            random.seed(0)
            random.shuffle(markers)
            
        # Cycle them in case there are not enough
        series_symbols = itertools.cycle(markers)
        
        # Work out the order to do the series in
        if options.categories is not None:
            # The user specified an order
            category_order = options.categories
        elif options.sort:
            # We need to sort the input categories ourselves
            category_order = sorted(series.iterkeys(), key=natural_keys)
        else:
            # Grab the series names in the order they originally appeared.
            category_order = initial_series_order.keys()
            
        if options.dotplot:
            # Assign X coordinates in series order
            for i, category in enumerate(category_order):
                # Each item gets an x coordinate equal to its series number
                series[category][0] = [i] * len(series[category][1])
    
        for series_name, series_color, series_symbol in \
            itertools.izip(category_order, series_colors, 
            series_symbols):
                
            # How do we want to plot (line or just scatter?)
            plot_func = pyplot.plot if options.lines else pyplot.scatter
                
            # Do the actual plot
            plot_func(series[series_name][0], series[series_name][1],
                label=series_name, color=series_color, marker=series_symbol)
                
    else:
        # Just plot the only series in the default color
        pyplot.scatter(series[""][0], series[""][1])
        
    # StackOverflow provides us with font sizing
    # <http://stackoverflow.com/q/3899980/402891>
    matplotlib.rcParams.update({"font.size": options.font_size})
    if options.show_n:
        # Add an n value to the title
        options.title += " (n = {})".format(sum((len(points)
            for (points, _) in series.itervalues())))
    pyplot.title(options.title)
    pyplot.xlabel(options.x_label)
    if options.log_x:
        # Turn on log X axis if desired. See
        # <http://stackoverflow.com/a/3513577/402891>
        pyplot.xscale("log")
    
    pyplot.ylabel(options.y_label)
    if options.log_y:
        # And log Y axis if desired.
        pyplot.yscale("log")
    
    if options.dotplot:
        # Turn off the x ticks
        pyplot.gca().get_xaxis().set_ticks([])
        # Set the plot bounds to just around the data
        pyplot.xlim((-1, len(category_order)))    
    
    # Apply any range restrictions
    if(options.min_x is not None):
        pyplot.xlim((options.min_x, pyplot.xlim()[1]))    
    if(options.max_x is not None):
        pyplot.xlim((pyplot.xlim()[0], options.max_x))
        
    if(options.min_y is not None):
        pyplot.ylim((options.min_y, pyplot.ylim()[1]))    
    if(options.max_y is not None):
        pyplot.ylim((pyplot.ylim()[0], options.max_y))
        
    if options.sparse_ticks:
        # Set up tickmarks to have only 2 per axis, at the ends
        pyplot.gca().xaxis.set_major_locator(
            matplotlib.ticker.FixedLocator(pyplot.xlim()))
        pyplot.gca().yaxis.set_major_locator(
            matplotlib.ticker.FixedLocator(pyplot.ylim()))
    
    # Make sure tick labels don't overlap. See
    # <http://stackoverflow.com/a/20599129/402891>
    pyplot.gca().tick_params(axis="x", pad=0.5 * options.font_size)
    
    # Make everything fit
    pyplot.tight_layout()
    
    if use_series:
        # Add a legend if we have multiple series
        
        if options.legend_overlay is None:
            # We want the default legend, off to the right of the plot.
        
            # First shrink the plot to make room for it.
            # TODO: automatically actually work out how big it will be.
            bounds = pyplot.gca().get_position()
            pyplot.gca().set_position([bounds.x0, bounds.y0, bounds.width * 0.5, 
                bounds.height])
                
            # Make the legend
            pyplot.legend(loc="center left", bbox_to_anchor=(1, 0.5))
            
        else:
            # We want the legend on top of the plot at the user-specified
            # location, and we want the plot to be full width.
            pyplot.legend(loc=options.legend_overlay)
    
    if options.save is not None:
        # Save the figure to a file
        pyplot.savefig(options.save, dpi=options.dpi)
    else:
        # Show the figure to the user
        pyplot.show()
        
    return 0

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
