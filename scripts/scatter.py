#!/usr/bin/env python2.7
"""
scatter: plot a scatterplot of a file of numbers. Numbers should be floats, two 
per line.

Re-uses sample code and documentation from 
<http://users.soe.ucsc.edu/~karplus/bme205/f12/Scaffold.html>
"""

import argparse, sys, os, itertools, math, collections, random, re
import matplotlib, matplotlib.ticker, matplotlib.cm, numpy
import copy

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
    parser.add_argument("--sparse_axes", action="store_true",
        help="use only bottom and left axes")
    parser.add_argument("--lines", action="store_true",
        help="connect points together")
    parser.add_argument("--tsv", action="store_true",
        help="use only tabs as separators in input file")
    parser.add_argument("--no_sort", dest="sort", action="store_false",
        help="do not sort categories from input file")
    parser.add_argument("--categories", nargs="+", default=None,
        help="categories to plot, in order")
    parser.add_argument("--category_labels", nargs="+",
        default=None,
        help="labels for all categories, in order")
    parser.add_argument("--legend_overlay", default=None,
        help="display the legend overlayed on the graph at this location")
    parser.add_argument("--no_legend", dest="show_legend", action="store_false",
        help="don't draw legend")
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
    parser.add_argument("--marker_size", type=float, default=None,
        help="marker size")
    parser.add_argument("--line_width", type=float, default=None,
        help="line_width")
    parser.add_argument("--annotate", action="store_true",
        help="annotate every point with its category name")
    
    return parser.parse_args(args)
    
def physics_layout_labels(to_label, series, other_spring=0.06,
    other_dist = 0.3, data_spring=0.02, data_dist = 0.2, target_spring=0.05,
    target_dist=0.15, max_steps=1000, min_x = 0, min_y = 0, max_x = 1,
    max_y = 1):
    """
    Given a series dict of points to label by series name, then a list for x or
    y, and then in a list by point number, and a similar structure of points to
    avoid, return a similar dict with the best locations for text to annotate
    those points.
    
    other_spring and other_dist control spring force for avoiding other points
    being laid out. data_spring and data_dist control spring force for avoiding
    data points. target_spring and target_dist control spring force for seeking
    each laid out point's own data point.
    
    max_steps controls how many iterations to run for.
    
    min_x, min_y, max_x, and max_y set a bounding box that points are forced to
    stay in.
    
    """
    
    # Deep copy the points to be labeled so we can update in place to move
    # things.
    positions = copy.deepcopy(to_label)
    
    # Keep forces that we accumulate into, organized the same way
    forces = collections.defaultdict(lambda: [[], []])
    
    # Set up bounds by dimension number
    min_bounds = [min_x, min_y]
    max_bounds = [max_x, max_y]
    
    # Compute x and y dimensions
    x_width = max_x - min_x
    y_height = max_y - min_y
    
    # We're going to do our layout in an imaginary 1 by 1 space within these
    # bounds, and convert vectors into it (divide by width or height) when doing
    # spring math and out of it (multiply by width or height) when applying
    # forces.
    
    def apply_forces():
        """
        Apply all the forces and make changes in positions. Also enforces the
        bounding box.
        
        """
        
        for name, dimensions in forces.iteritems():
            for dimension, values in enumerate(dimensions):
                for i in xrange(len(values)):
                    # Use each force as an offset
                    
                    # First make a new number that we'll update
                    new_pos = positions[name][dimension][i]
                    
                    # Update it
                    new_pos += values[i]
                    
                    # Don't let it escape the box
                    if new_pos < min_bounds[dimension]:
                        new_pos = min_bounds[dimension]
                    if new_pos > max_bounds[dimension]:
                        new_pos = max_bounds[dimension]
                        
                    # Apply it
                    positions[name][dimension][i] = new_pos

    def reset_forces():
        """
        Clear all forces to 0 in both dimensions
        """
        for name, dimensions in positions.iteritems():
            for dimension, values in enumerate(dimensions):
                forces[name][dimension] = [0.0] * len(values)
                
    def apply_force(name, index, x, y):
        """
        Apply the given force in x and y to the given point number in the given
        series.
        """
        
        forces[name][0][index] += x
        forces[name][1][index] += y
                
    def for_each_point(points_dict):
        """
        Go through each point in a positions or forces dict and yield (series,
        index, x value, y value).
        
        """
        
        for name, dimensions in points_dict.iteritems():
            for i in xrange(len(dimensions[0])):
                yield (name, i, dimensions[0][i], dimensions[1][i])
        
    
    # Start with 0 forces
    reset_forces()
    
    for step in xrange(max_steps):
        # Do a bunch of steps
    
        for point_series, point_index, point_x, point_y in \
            for_each_point(positions):        
            # For each text point
            
            for other_series, other_index, other_x, other_y in \
                for_each_point(positions):
                # Get its offset from each other text point and apply an away
                # force when too close
                
                if point_series == other_series and point_index == other_index:
                    # Don't affect self
                    continue
                
                # What's the offset in each dimension, in spring space?
                x_offset = (other_x - point_x) / x_width
                y_offset = (other_y - point_y) / y_height
                
                # And the total offset
                offset_length = math.pow(math.pow(x_offset, 2) +
                    math.pow(y_offset, 2), 0.5)
                
                if offset_length < other_dist:
                    # Too close! Spring away!
                    diff = other_dist - offset_length
                    
                    if offset_length == 0:
                        # Go in a random not-here direction
                        # TODO: biased towards corners
                        x_offset = random.random() - 0.5
                        if x_offset == 0:
                            # If we actually hit 0, move
                            x_offset = 0.1
                        y_offset = random.random() - 0.5
                        offset_length = math.pow(math.pow(x_offset, 2) +
                            math.pow(y_offset, 2), 0.5)
                    
                    x_force = -x_offset / offset_length * other_spring * diff
                    y_force = -y_offset / offset_length * other_spring * diff
                    
                    apply_force(point_series, point_index, x_force * x_width,
                        y_force * y_height)
                    
            for data_series, data_index, data_x, data_y in \
                for_each_point(series):
                # Get its offset from each data point and apply an away
                # force when too close
                
                # What's the offset in each dimension
                x_offset = (data_x - point_x) / x_width
                y_offset = (data_y - point_y) / y_height
                
                # And the total offset
                offset_length = math.pow(math.pow(x_offset, 2) +
                    math.pow(y_offset, 2), 0.5)
                
                if offset_length < data_dist:
                    # Too close! Spring away!
                    diff = data_dist - offset_length
                    
                    if offset_length == 0:
                        # Go in a random not-here direction
                        # TODO: biased towards corners
                        x_offset = random.random()
                        if x_offset == 0:
                            # If we actually hit 0, move
                            x_offset = 0.1
                        y_offset = random.random()
                        offset_length = math.pow(math.pow(x_offset, 2) +
                            math.pow(y_offset, 2), 0.5)
                    
                    x_force = -x_offset / offset_length * data_spring * diff
                    y_force = -y_offset / offset_length * data_spring * diff
                    
                    apply_force(point_series, point_index, x_force * x_width,
                        y_force * y_height)
        
            # Get its offset from its data point and apply a toward force
            target_x = to_label[point_series][0][point_index]
            target_y = to_label[point_series][1][point_index]
            
            # What's the offset in each dimension
            x_offset = (target_x - point_x) / x_width
            y_offset = (target_y - point_y) / y_height
            
            # And the total offset
            offset_length = math.pow(math.pow(x_offset, 2) +
                math.pow(y_offset, 2), 0.5)
            
            if offset_length > target_dist:
                # Too far! Spring towards!
                diff = offset_length - target_dist
                
                x_force = x_offset / offset_length * target_spring * diff
                y_force = y_offset / offset_length * target_spring * diff
                
                apply_force(point_series, point_index, x_force * x_width,
                    y_force * y_height)
        
        # Apply all the forces
        apply_forces()
        
        # Zero forces
        reset_forces()
    
    # Return the updated positions
    return positions
    

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
        elif len(parts) >= 3:
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
            # Don't sort, use the input order
            options.sort = False
        
            # Fix up the options so we don't ask for categories with no points.
            for i in xrange(len(options.categories)):
                if series.has_key(options.categories[i]):
                    # Leave this one
                    continue
                    
                # Otherwise none it out from everything
                options.categories[i] = None
                
                if options.colors is not None:
                    options.colors[i] = None
                
                if options.markers is not None and len(options.markers) > i:
                    options.markers[i] = None
                    
                if options.category_labels is not None:
                    options.category_labels[i] = None
             
            # Now filter out everrything we noned out       
            options.categories = [x for x in options.categories
                if x is not None]
                
            if options.colors is not None:
                options.colors = [x for x in options.colors
                if x is not None]
                
            if options.markers is not None:
                options.markers = [x for x in options.markers
                if x is not None]
                
            if options.category_labels is not None:
                options.category_labels = [x for x in options.category_labels
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
            
        # Assign names to categories
        if options.category_labels is None:
            category_names = {category: category for category in category_order}
        else:
            category_names = {category: label for (category, label) in 
                itertools.izip(category_order, options.category_labels)}
            
        if options.dotplot:
            # Assign X coordinates in series order
            for i, category in enumerate(category_order):
                # Each item gets an x coordinate equal to its series number
                series[category][0] = [i] * len(series[category][1])
    
        print category_names
    
        for series_name, series_color, series_symbol in \
            itertools.izip(category_order, series_colors, 
            series_symbols):
                
            # How do we want to plot (line or just scatter?)
            plot_func = pyplot.plot if options.lines else pyplot.scatter

            # Build up some different options depending on plot style
            plot_opt = {}
            if options.lines:
                # black marker outline doesn't look great on lines
                plot_opt["markeredgecolor"] = series_color
                # line and marker size options:
                if options.marker_size is not None:
                    plot_opt["ms"] = options.marker_size
                if options.line_width is not None:
                    plot_opt["linewidth"] = options.line_width
            elif options.marker_size is not None:
                plot_opt["s"] = options.marker_size
                    
            # Do the actual plot
            plot_func(series[series_name][0], series[series_name][1],
                label=category_names[series_name], color=series_color,
                      marker=series_symbol, **plot_opt)
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
        
    if use_series and options.annotate:
        # We need to annotate every point in every series with its
        # series name, trying not to overlap points or other annotations.
        
        # Figure out where to put the labels. We want them clear of the axes.
        # How wide is the plot?
        x_width = pyplot.xlim()[1] - pyplot.xlim()[0]
        y_height = pyplot.ylim()[1] - pyplot.ylim()[0]
        
        # Retract a certain distance form the edges
        min_x = pyplot.xlim()[0] + 0.1 * x_width
        max_x = pyplot.xlim()[1] - 0.1 * x_width
        min_y = pyplot.ylim()[0] + 0.1 * y_height
        max_y = pyplot.ylim()[1] - 0.1 * y_height
        
        # Find the center point of each series and decide to label it
        to_label = {}
        for series_name, series_points in series.iteritems():
            # For each series, how many points are in it?
            series_length = len(series_points[0])
            # What's the index of the middle point?
            series_center = int(series_length/2)
            
            # Clip the series down to the central point and copy it over.
            to_label[series_name] = [[series_points[0][series_center]],
                [series_points[1][series_center]]]
        
        # Pass the whole dict from name to list of x/y lists into the physics
        # layout within that box. We'll get a similar dict back with the best
        # label position for each point.
        label_positions = physics_layout_labels(to_label, series,
            min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
        
        # Reset series color iterator to the beginning
        series_colors = itertools.cycle(colors)
        
        for series_name, series_color in itertools.izip(category_order,
            series_colors):
            
            # For each series
            for i in xrange(len(label_positions[series_name][0])):
                # For each label position
        
                # Label the point with an arrow in the correct color
                # Make sure to center the text on the text position
                pyplot.gca().annotate(category_names[series_name], 
                    xy=(to_label[series_name][0][i],
                    to_label[series_name][1][i]), 
                    xytext=(label_positions[series_name][0][i],
                    label_positions[series_name][1][i]),
                    color=series_color,
                    horizontalalignment="center",
                    verticalalignment="center",
                    arrowprops=dict(arrowstyle="->", color=series_color))
        
        
    if options.sparse_ticks:
        # Set up tickmarks to have only 2 per axis, at the ends
        pyplot.gca().xaxis.set_major_locator(
            matplotlib.ticker.FixedLocator(pyplot.xlim()))
        pyplot.gca().yaxis.set_major_locator(
            matplotlib.ticker.FixedLocator(pyplot.ylim()))
            
    if options.sparse_axes:
        # Don't draw top or right axes
        pyplot.gca().spines["right"].set_visible(False)
        pyplot.gca().spines["top"].set_visible(False)
        # Or their tick marks
        pyplot.gca().yaxis.set_ticks_position("left")
        pyplot.gca().xaxis.set_ticks_position("bottom")
    
    # Make sure tick labels don't overlap. See
    # <http://stackoverflow.com/a/20599129/402891>
    pyplot.gca().tick_params(axis="x", pad=0.5 * options.font_size)
    
    # Make everything fit
    pyplot.tight_layout()
    
    if use_series and options.show_legend:
        # Add a legend if we have multiple series
        
        if options.legend_overlay is None:
            # We want the default legend, off to the right of the plot.
        
            # First shrink the plot to make room for it.
            # TODO: automatically actually work out how big it will be.
            bounds = pyplot.gca().get_position()
            pyplot.gca().set_position([bounds.x0, bounds.y0, bounds.width * 0.5, 
                bounds.height])
                
            # Make the legend
            pyplot.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))
            
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
