#!/usr/bin/env python2.7
"""
barchart: plot a bar chart of a TSV file of numbers. The file should be a
column of int or text labels and a column of floats, with one value per label.

Re-uses sample code and documentation from 
<http://users.soe.ucsc.edu/~karplus/bme205/f12/Scaffold.html>
"""

import argparse, sys, os, itertools, math, collections, re
import matplotlib, matplotlib.ticker

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
    parser.add_argument("--divisions", type=argparse.FileType('r'),
        default=[],
        help="file of numbers to be counted towards a lower bar instead")
    parser.add_argument("--title", default="Bar Chart",
        help="the plot title")
    parser.add_argument("--x_label", default="X",
        help="the X aximatplotlib.tickers label")
    parser.add_argument("--y_label", default="Y",
        help="the Y axis label")
    parser.add_argument("--log_y", action="store_true",
        help="log Y axis")
    parser.add_argument("--font_size", type=int, default=12,
        help="the font size for text")
    parser.add_argument("--save",
        help="save figure to the given filename instead of showing it")
    parser.add_argument("--dpi", type=int, default=300,
        help="save the figure with the specified DPI, if applicable")
    parser.add_argument("--sparse_ticks", action="store_true",
        help="Use sparse tick marks")
    parser.add_argument("--x_sideways", action="store_true",
        help="write X axis labels vertically")
    parser.add_argument("--bar_width", type=float, default=0.8,
        help="width of bars as a fraction of maximum")
    parser.add_argument("--min", type=float, default=None,
        help="minimum value allowed")
    parser.add_argument("--max", type=float, default=None,
        help="maximum value allowed")
    parser.add_argument("--categories", nargs="+", default=None,
        help="categories to plot, in order")
    parser.add_argument("--category_labels", nargs="+", default=None,
        help="labels for all categories, in order")
    parser.add_argument("--colors", nargs="+", default=None,
        help="colors for all categories, in order")
        
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
    
    # Read the divisions, if applicable, and store them by category name. We
    # will take so much from the normal bar and allocate it to a bar that
    # appears below it.
    divisions = collections.defaultdict(lambda: 0)
    for line in options.divisions:
        # Unpack and parse the two numbers on this line (category and value)
        parts = line.strip().split('\t')
        if len(parts) < 2:
            # Skip empty/invalid lines
            continue
            
        try:
            # Parse categories to ints if possible
            category = int(parts[0])
        except ValueError:
            category = parts[0]
        value = float(parts[1])
        
        # Sum in the values
        divisions[category] += value
    
    
    # This dict holds the value for every bar. Each starts at 0.
    categories = collections.defaultdict(lambda: 0)
    
    for line in options.data:
        # Unpack and parse the two numbers on this line (category and value)
        parts = line.strip().split('\t')
        if len(parts) < 2:
            # Skip empty/invalid lines
            continue
            
        try:
            # Parse categories to ints if possible
            category = int(parts[0])
        except ValueError:
            category = parts[0]
        value = float(parts[1])
        
        # Sum in the values
        categories[category] += value
        
    for category, division in divisions.iteritems():
        # Allocate space for the lower bar.
        categories[category] -= division
    
    # This holds the category values after min/max filtering.
    categories_filtered = {}
    
    for category, value in categories.iteritems():
        
        if options.min is not None and value < options.min:
            # Throw out values that are too small (like 0s for log_y)
            continue
            
        if options.max is not None and value > options.max:
            # Throw out values that are too large
            continue
            
        # If the total value in a category passes the test, keep it
        categories_filtered[category] = value
        
    # Throw away the original data and keep the filtered data.
    categories = categories_filtered
    
    # Compute the order to use
    if options.categories is not None:
        category_order = options.categories
    else:
        category_order = sorted(categories.iterkeys(), key=natural_keys)
    
    # Work out what category labels to use
    category_labels = options.category_labels \
        if options.category_labels is not None else category_order
        
    # Work out what colors to use
    category_colors = options.colors \
        if options.colors is not None else ['b'] * len(category_order)
   
    for i in xrange(len(category_order)):
        # Look at every category
        if not categories.has_key(category_order[i]):
            # We need to remove this one since we have no data for it
            # We do this by none-ing it out in all the lists and then filtering
            category_order[i] = None
            category_labels[i] = None
            category_colors[i] = None
            
    # Do the filtering of things we noned out.
    category_order = [x for x in category_order if x is not None]
    category_labels = [x for x in category_labels if x is not None]
    category_colors = [x for x in category_colors if x is not None]
    
    # Do the plot
    pyplot.bar(xrange(len(category_order)), [categories[category] for category in 
        category_order], bottom=[divisions[category] for category in 
        category_order], color=category_colors, width=options.bar_width)
    # Plot the below-division bars.
    pyplot.bar(xrange(len(category_order)), [divisions[category] for category in 
        category_order], width=options.bar_width)
        
    # StackOverflow provides us with font sizing
    # http://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
    matplotlib.rcParams.update({"font.size": options.font_size})
    pyplot.title("{} (n = {})".format(options.title, len(category_order)))
        
    pyplot.xlabel(options.x_label)
    
    # Label the columns with the appropriate text. Account for 1-based ticks.
    pyplot.xticks([x + (options.bar_width/2) for x in xrange(len(categories.keys()))],
        category_labels, rotation=90 if options.x_sideways else 0)
    
    pyplot.ylabel(options.y_label)
    if options.log_y:
        # And log Y axis if desired.
        pyplot.yscale("log")
        
    if options.max is not None:
        # Set only the upper y limit
        pyplot.ylim((pyplot.ylim()[0], options.max))
        
    if options.sparse_ticks:
        # Set up tickmarks to have only 2 per axis, at the ends
        pyplot.gca().yaxis.set_major_locator(
            matplotlib.ticker.FixedLocator(pyplot.ylim()))
            
    # Get rid of empty space on the right, by moving the right end of the plot
    # to the right end of the last bar.
    pyplot.xlim(0, len(category_order) - 1 + options.bar_width)
    
    # Make sure tick labels don't overlap. See
    # <http://stackoverflow.com/a/20599129/402891>
    pyplot.gca().tick_params(axis="x", pad=0.5 * options.font_size)
    
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
