#!/usr/bin/env python2.7
"""
boxplot: plot a boxplot of a TSV file of numbers. The file should be a
column of int or text labels and a column of floats.

Re-uses sample code and documentation from 
<http://users.soe.ucsc.edu/~karplus/bme205/f12/Scaffold.html>
"""

import argparse, sys, os, itertools, math, collections, re
import matplotlib, matplotlib.ticker, matplotlib.lines, numpy

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


def flatten(iterables):
    """
    Flatten an iterable of iterables into a single iterable. See
    <https://docs.python.org/2/library/itertools.html>.
    
    """
    return itertools.chain.from_iterable(iterables)

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
    parser.add_argument("--title", default="Box Plot",
        help="the plot title")
    parser.add_argument("--x_label", default="X",
        help="the X axis label")
    parser.add_argument("--y_label", default="Y",
        help="the Y axis label")
    parser.add_argument("--log_y", action="store_true",
        help="log Y axis")
    parser.add_argument("--hline", type=float, default=None,
        help="draw a horizontal line at the given Y value")
    parser.add_argument("--hline_median", default=None,
        help="draw a horizontal line at the median of the given category")
    parser.add_argument("--font_size", type=int, default=12,
        help="the font size for text")
    parser.add_argument("--save",
        help="save figure to the given filename instead of showing it")
    parser.add_argument("--dpi", type=int, default=300,
        help="save the figure with the specified DPI, if applicable")
    parser.add_argument("--line_width", type=float, default=1,
        help="the width of lines making up the boxes and their parts")
    parser.add_argument("--sparse_ticks", action="store_true",
        help="use sparse tick marks")
    parser.add_argument("--x_sideways", action="store_true",
        help="write X axis labels vertically")
    parser.add_argument("--min", type=float, default=None,
        help="minimum Y allowed")
    parser.add_argument("--max", type=float, default=None,
        help="maximum Y value")
    parser.add_argument("--max_max", type=float, default=None,
        help="limit on maximum Y value")
    parser.add_argument("--means", action="store_true",
        help="include means for each category")
    parser.add_argument("--no_n", dest="show_n", action="store_false",
        help="don't add n value to title")
        
    # We take these next 4 options, optionally, in strided groupings. So you can
    # have multiple --categories options, each of which is named with a
    # --grouping and has its categories labeled with a --category_labels
    parser.add_argument("--categories", nargs="+", action="append",
        default=None,
        help="categories to plot, in order")
    parser.add_argument("--category_labels", nargs="+", action="append",
        default=None,
        help="labels for all categories, in order")
    parser.add_argument("--colors", nargs="+", action="append",
        default=None,
        help="colors for all categories, in order")
    parser.add_argument("--grouping", action="append", default=None,
        help="name for corresponding --categories/--category_labels options")
    
    parser.add_argument("--grouping_colors", nargs="+", default=None,
        help="Matplotlib colors for all the groupings")
    
    # If we have colors and names for groupings, we will use a legend.
    parser.add_argument("--legend_overlay", default=None,
        help="display any legend overlayed on the graph at this location")
    parser.add_argument("--legend_columns", type=int, default=1,
        help="number of columns to use in any legend")
    parser.add_argument("--no_legend", dest="show_legend", action="store_false",
        help="don't add a legend")
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
    
    # Make the figure with the appropriate size
    pyplot.figure(figsize=(options.width, options.height))
    
    # This defaultdict holds one list for every category, containing all the
    # floats in that category.
    categories = collections.defaultdict(list)
    
    max_found = None
    min_found = None
    
    for line in options.data:
        # Unpack and parse the two numbers on this line (category and value)
        parts = line.strip().split('\t')
        if len(parts) < 2:
            # Skip empty/invalid lines
            continue
         
        category = parts[0]
        
        try:
            # Parse categories to ints if possible, for numerical sort
            category = int(category)
        except ValueError:
            pass
            
        value = float(parts[1])
        
        if options.min is not None and value < options.min:
            # Throw out values that are too small (like 0s for log_y)
            continue
            
        if options.max is not None and value > options.max:
            # Throw out values that are too large
            continue
            
        if max_found is None or value > max_found:
            # Track the max value we find
            max_found = value
            
        if min_found is None or value < min_found:
            # Track the min value we find
            min_found = value
        
        # Put each in the appropriate list.
        categories[category].append(value)
    
    
    
    # Should we use groupings?
    if ((options.categories is not None and len(options.categories) > 0) or 
        (options.category_labels is not None and
        len(options.category_labels) > 0) or 
        (options.colors is not None and len(options.colors) > 0)):
        
        # The user specified some option that implies we are supposed to use
        # groupings. Don't check consistency, but do go into groupings mode.
        use_groupings = True
    else:
        # Don't use groupings
        use_groupings = False
        
    # Compute the order to use
    if options.categories is not None:
        # We were given a list of lists of categories. Flatten them.
        category_order = list(flatten(options.categories))
    else:
        # Make our own order by natural sort.
        category_order = sorted(categories.iterkeys(), key=natural_keys)
    
    # Work out category labels
    if options.category_labels is not None:
        # The user specified labels
        category_labels = list(flatten(options.category_labels))
        
        if len(category_labels) != len(category_order):
            print(category_labels)
            print(category_order)
            # Complain that the user doesn't know what they want.
            raise Exception(
                "Incorrect number of labels ({}) for number of "
                "categories ({})".format(len(category_labels),
                len(category_order)))
    else:
        # Just label each category with its actual name.
        category_labels = category_order
        
    if options.colors:
        # Make a single list of colors
        category_colors = list(flatten(options.colors))
    else:
        # Automatic colors
        category_colors = None
        
    for i in xrange(len(category_order)):
        # Look at every category
        if not categories.has_key(category_order[i]):
            # We need to remove this one since we have no data for it
            # We do this by none-ing it out in all the lists and then filtering
            category_order[i] = None
            category_labels[i] = None
            if category_colors is not None:
                category_colors[i] = None
            
    # Do the filtering of things we noned out.
    category_order = [x for x in category_order if x is not None]
    category_labels = [x for x in category_labels if x is not None]
    if category_colors is not None:
        category_colors = [x for x in category_colors if x is not None]
    
    if options.means:
        # We want to see mean values.
        for i in xrange(len(category_order)):
            # Calculate the mean of each category and add it to the label.
            category_labels[i] += "\n(mean={:.3f})".format(sum(
                categories[category_order[i]]) / float(len(
                categories[category_order[i]])))
                
    # We have a dict that we add all our per-column properties to. TODO: needs
    # new-ish matplotlib
    base_props = {'linewidth': options.line_width}
        
    if use_groupings:
        # We have groupings to group.
        
        # What are the starts and lengths of each grouping in the list of
        # categories? Store as a (start, length) tuple list.
        grouping_bounds = []
        
        # Where does the next grouping start?
        next_start = 0
        
        for grouping_members in (options.categories
            if options.categories is not None else options.category_labels):
            # For each list of categories in a grouping, pulling from any option
            # that reveals grouping structure...
            
            # Throw out the members that don't exist
            grouping_members = [m for m in grouping_members
                if categories.has_key(m)]
            
            # Work out the start category and length, and record it
            grouping_bounds.append((next_start, len(grouping_members)))
            
            # Next grouping starts where this one ended
            next_start += len(grouping_members)
            
        for i, (start, length) in enumerate(grouping_bounds):
            # For each grouping with its bounds in the list of categories...
            
            if options.grouping is not None and (options.colors is not None or
                options.grouping_colors is None):
                # Stick the grouping name on the first category label as a
                # second line, since we can't use a legend.
                category_labels[start] += "\n" + options.grouping[i]

            # Pull out the values that need to be plotted
            values = [categories[category] for category in 
                category_order[start:start + length]]
                
            # And figure out where the bars go.
            positions = [x + 1 for x in xrange(start, start + length)]

            for j in xrange(length):
                # For each box we have to put in the grouping
    
                # We need to set its properties
                boxprops = dict(base_props)
                
                if options.colors is not None:
                    # We need to do individual item colors
                    boxprops["color"] = category_colors[start + j]
                    
                    # Do the plot for this single box, at the correct place,
                    # passing the color options if needed, as well as all the
                    # other options.
                    pyplot.boxplot([values[j]], positions=[positions[j]],
                        widths=0.5, boxprops=boxprops, meanprops=boxprops,
                        medianprops=boxprops, flierprops=boxprops,
                        whiskerprops=boxprops, capprops=boxprops)
                    
                elif options.grouping_colors is not None:
                    # Set a color for this whole grouping
                    boxprops["color"] = options.grouping_colors[i]
                    
                    # Only color the box itself to match old behavior.
                    pyplot.boxplot([values[j]], positions=[positions[j]],
                        widths=0.5, boxprops=boxprops, meanprops=base_props,
                        medianprops=base_props, flierprops=base_props,
                        whiskerprops=base_props, capprops=base_props)
                    
                
                
        # We need to fix up the x limits since boxplot messes them up.
        pyplot.xlim(0.5, len(category_labels) + 0.5)
                
    else:
        # Do the plot for all the categories at once.
        
        pyplot.boxplot([categories[category] for category in category_order],
            boxprops=base_props, meanprops=base_props, medianprops=base_props,
            whiskerprops=base_props, capprops=base_props)
            
            
        
    if options.hline is not None:
        # Add in our horizontal line that the user asked for.
        pyplot.axhline(y=options.hline, color='r', linestyle='--')
        
    if options.hline_median is not None:
        # Add in our horizontal mean line
        
        # Compute the mean
        category_median = numpy.median(categories[options.hline_median])
        
        # Find the color
        line_color = 'r'
        for i in xrange(len(category_colors)):
            if category_order[i] == options.hline_median:
                line_color = category_colors[i]
                break
        
        pyplot.axhline(y=category_median, color=line_color, linestyle='--')
        
        # We want to mark the best thing with its percent deviation
        best_category = None
        best_deviation = None
        
        for i, category in enumerate(category_order):
            # Apply a +/- median difference percentage to the label
            
            other_median = numpy.median(categories[category])
            
            portion = other_median / category_median
            
            percent = (portion - 1) * 100
            
            if best_deviation is None or percent > best_deviation:
                # We found the best thing
                best_category = i
                best_deviation = percent
        
        if best_category is not None:
            # Apply a +/-% label to the best thing
            category_labels[best_category] += "\n({:+.2f}%)".format(
                best_deviation)
        
    # StackOverflow provides us with font sizing
    # http://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
    matplotlib.rcParams.update({"font.size": options.font_size})
    if options.show_n:
        # Add an n value to the title
        options.title += " (n = {})".format(sum((len(categories[category])
            for category in category_order)))
    pyplot.title(options.title)
        
        
    pyplot.xlabel(options.x_label)
    # Label the columns with the appropriate text. Account for 1-based ticks.
    pyplot.xticks(xrange(1, len(category_labels) + 1), category_labels, 
        rotation=90 if options.x_sideways else 0)
    
    if options.max_max < max_found:
        # Bring in the upper limit
        pyplot.ylim((pyplot.ylim()[0], options.max_max))
    
    if options.min is not None:
        # Set only the lower y limit
        pyplot.ylim((options.min, pyplot.ylim()[1]))
    if options.max is not None:
        # Set only the upper y limit
        pyplot.ylim((pyplot.ylim()[0], options.max))
    
    pyplot.ylabel(options.y_label)
    if options.log_y:
        # And log Y axis if desired.
        pyplot.yscale("log")
        
    if options.sparse_ticks:
        # Set up tickmarks to have only 2 per axis, at the ends
        pyplot.gca().yaxis.set_major_locator(
            matplotlib.ticker.FixedLocator(pyplot.ylim()))
    
    # Make sure tick labels don't overlap. See
    # <http://stackoverflow.com/a/20599129/402891>
    pyplot.gca().tick_params(axis="x", pad=0.5 * options.font_size)
    
    # Make everything fit
    pyplot.tight_layout()
    
    
    if (use_groupings and options.grouping_colors is not None and
        options.grouping is not None) and options.show_legend:
        # We have groupings, and we can make a legend. We need to do it here
        # since we need to resize the plot after it is tightly laid out.
        
        # We're going to make proxy artists and a legend if our groupings
        # have names and colors. See
        # <http://matplotlib.org/users/legend_guide.html#proxy-legend-
        # handles>
        
        # Make a labeled artist for each grouping
        artists = [matplotlib.lines.Line2D([], [],
            color=options.grouping_colors[i], label=options.grouping[i]) 
            for i in xrange(len(options.grouping))]
        
        # Make a dict to configure the legend
        legend_config = {"handles": artists, "ncol": options.legend_columns,
            "loc": options.legend_overlay}
    
    
        if options.legend_overlay is None:
            # We want the default legend, off to the right of the plot.
        
            # First shrink the plot to make room for it.
            # TODO: automatically actually work out how big it will be.
            bounds = pyplot.gca().get_position()
            pyplot.gca().set_position([bounds.x0, bounds.y0,
                bounds.width * 0.5, bounds.height])
                
            # Make the legend be on the side
            legend_config["loc"] = "center left"
            legend_config["bbox_to_anchor"] = (1, 0.5)
        
    
        # Make the legend we configured, passing all those keyword arguments
        pyplot.legend(**legend_config)
    
    if options.save is not None:
        # Save the figure to a file
        pyplot.savefig(options.save, dpi=options.dpi)
    else:
        # Show the figure to the user
        pyplot.show()
        
    return 0

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
