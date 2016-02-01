__author__ = 'Jwely'

import os
from py.tex.TeXWriter import TeXWriter
from py.piv import shorthand_to_tex
from py.piv.construct_experiments import construct_experiments
from py.config import *


class TeXRunFigurePage(TeXWriter):
    def __init__(self, main_path, tex_title, experiment_id, force_recalc=False):
        """
        This class uses constructors to build up a data set from the raw v3d files, perform all processing,
        create figures, and compile a publication ready .tex file from all the arguments. It places all outputs
        in the texdocs/figs directory. All tex files in the texdocs/figs/tex have their corresponding figures

        :param main_path:           filepath to main.tex
        :param tex_title:           title for this tex file. will be placed in texdocs/figs/tex
        :param experiment_id:       an integer experiment ID from which to synthesize plots
        :param force_recalc:
        :return:
        """

        # invoke parent class init and add a few other attributes
        main_dir = os.path.dirname(os.path.abspath(main_path))
        self.experiment_id = experiment_id
        self.tex_title = tex_title
        self.figure_dir = os.path.join(main_dir, "figs/{0}".format(self.tex_title))
        if not os.path.exists(self.figure_dir):
            os.mkdir(self.figure_dir)
        TeXWriter.__init__(self, main_path, os.path.join(main_dir, "figs/tex/{0}.tex".format(self.tex_title)))

        # ensure inputs are properly formatted
        if not isinstance(experiment_id, int):
            raise Exception("please input experiment ID as an integer")

        self.experiment = construct_experiments(experiment_table_path=EXPERIMENT_TABLE_PATH,
                                                experiment_directory_path=DATA_FULL_DIR,
                                                ids=[experiment_id],
                                                min_points=DEFAULT_MIN_POINTS,
                                                force_recalc=force_recalc)[0]
        self.axial_vortex = self.experiment.axial_vortex


    def add_contour_plot(self, component, caption, width, special_tag=None, create_kwargs=None):
        """
        wraps TeXWriter.add_figure and the AxialVortex.contour_plot functions to add a contour
        plot figure to the TeXRunFigurePage with all the bells and whistles taken care of
        automatically. This is used to generate all the contour plots of interest for a given
        experiment run.

        :param component:       shorthand component to add contour plot
        :param caption:         Tex caption to add to figure
        :param width:           the width of the plot on the page such as '5in'
        :param create_kwargs:   manual kwargs for the contour_plot function
        :return:
        """

        if create_kwargs is None:
            create_kwargs = {}
        if special_tag is None:
            special_tag = ""
        if 'title' not in create_kwargs.keys():
            create_kwargs["title"] = shorthand_to_tex(component)
        create_kwargs["component"] = component

        create_from_function = self.axial_vortex.contour_plot

        figure_filename = "{0}_{1}_{2}_contour.png".format(self.tex_title, component, special_tag)
        figure_path = os.path.join(self.figure_dir, figure_filename)
        TeXWriter.add_figure(self, figure_path, caption, width, create_from_function, create_kwargs)


    def add_scatter_plot(self, component_x, component_y, caption, width,
                         special_tag=None, create_kwargs=None):
        """
        wraps TeXWriter.add_figure and the AxialVortex.scatter_plot functions to add a contour
        plot figure to the TeXRunFigurePage with all the bells and whistles taken care of
        automatically. This is used to generate all the contour plots of interest for a given
        experiment run.

        :param component_x:         x axis component to plot
        :param component_y:         y axis component to plot
        :param caption:             Tex caption to add to figure
        :param width:               the width of the plot on the page such as '5in'
        :param special_tag:         extra tag to place inside the filename for unique id of something with
                                    unique kwargs.
        :param create_kwargs:       manual kwargs for the scatter_plot function
        :return:
        """

        if create_kwargs is None:
            create_kwargs = {}
        if special_tag is None:
            special_tag = ""
        if 'title' not in create_kwargs.keys():
            create_kwargs["title"] = "{0} vs {1}".format(
                shorthand_to_tex(component_y), shorthand_to_tex(component_x))
        create_kwargs['component_x'] = component_x
        create_kwargs['component_y'] = component_y

        create_from_function = self.axial_vortex.scatter_plot

        figure_filename = "{0}_{1}vs{2}_{3}_scatter.png".format(
            self.tex_title, component_y, component_x, special_tag)

        figure_path = os.path.join(self.figure_dir, figure_filename)
        TeXWriter.add_figure(self, figure_path, caption, width, create_from_function, create_kwargs)


    def add_stream_plot(self, caption, width):
        """
        Adds a stream plot to the texdoc, has no special arguments.

        :param caption:         Tex caption to add to figure
        :param width:           the width of the plot on the page such as '5in'
        """

        create_kwargs = {}  # no interesting kwargs at the moment
        create_from_function = self.axial_vortex.stream_plot

        figure_filename = "{0}_stream.png".format(self.tex_title)
        figure_path = os.path.join(self.figure_dir, figure_filename)
        TeXWriter.add_figure(self, figure_path, caption, width, create_from_function, create_kwargs)


    def add_quiver_plot(self, caption, width):
        """
        Adds a quiver plot to the texdoc, has no special arguments.

        :param caption:         Tex caption to add to figure
        :param width:           the width of the plot on the page such as '5in'
        """

        create_kwargs = {}  # no interesting kwargs at the moment
        create_from_function = self.axial_vortex.quiver_plot

        figure_filename = "{0}_quiver.png".format(self.tex_title)
        figure_path = os.path.join(self.figure_dir, figure_filename)
        TeXWriter.add_figure(self, figure_path, caption, width, create_from_function, create_kwargs)


    def add_dynamic_plot(self, component_y, caption, width, special_tag=None, create_kwargs=None):
        """
        wraps TeXWriter.add_figure and the AxialVortex.dynamic functions to add a dynamic plot
        to the page

        :param component_y:         y axis component to plot
        :param caption:             Tex caption to add to figure
        :param width:               the width of the plot on the page such as '5in'
        :param special_tag:         extra tag to place inside the filename for unique id of something with
                                    unique kwargs.
        :param create_kwargs:       manual kwargs for the scatter_plot function
        :return:
        """

        if create_kwargs is None:
            create_kwargs = {}
        if special_tag is None:
            special_tag = ""
        if 'title' not in create_kwargs.keys():
            create_kwargs["title"] = "{0} over time".format(shorthand_to_tex(component_y))
        create_kwargs['component_y'] = component_y

        create_from_function = self.axial_vortex.dynamic_plot

        figure_filename = "{0}_{1}_{2}_dynamic.png".format(self.tex_title, component_y, special_tag)
        figure_path = os.path.join(self.figure_dir, figure_filename)
        TeXWriter.add_figure(self, figure_path, caption, width, create_from_function, create_kwargs)
