__author__ = 'Jwely'

import os
import traceback


class TeXFigureGenerator():


    def __init__(self, main_path, figure_filepath, caption, width):
        """
        Creates a list of Tex formatted strings to display the input relative filepath as
        a figure. Can call method `get_tex` to retrieve the list of strings, which can then be writen
        to an actual file line by line

        :param main_path:          filepath the tex document which is referencing this figure
        :param figure_filepath:    filepath to the image
        :param caption:            caption to use on the figure
        :param width:              width of the figure as a string, such as '5in'
        """

        # remove the extension for TeX \includegraphics command
        self.main_path = main_path
        self.absolute_filepath = os.path.abspath(figure_filepath)
        self.relative_filepath = self._path_relative_to_main(main_path, figure_filepath)

        # remove extension and isolate the filename to use as figure label
        self.label = "fig:{0}".format("".join(os.path.basename(self.relative_filepath).split(".")[:-1]))

        def enf_string(arg):

            if isinstance(arg, list) and len(arg) == 1:
                return str(arg[0])
            else:
                return str(arg)

        self.caption = enf_string(caption)
        self.width = enf_string(width)

    @staticmethod
    def _path_relative_to_main(main_path, figure_filepath):
        """
        converts a filepath to a path relative to the main_path, the location of main.tex

        :param figure_filepath: a filepath, absolute or relative to this script
        :return: the input filepath, but relative to the main_path
        """

        filepath = os.path.abspath(figure_filepath)
        if not os.path.exists(main_path):
            Warning("filepath '{0}' does not exist".format(filepath))
        common_prefix = os.path.commonprefix([filepath, main_path])
        return os.path.relpath(filepath, common_prefix)


    def create_from(self, function, kwargs=None):
        """
        Actually creates the figure from some input class instance and method and its arguments.
        Just wraps the figure creation arguments to ensure absolute_filepath points to a real place.
        input function must have "outpath" argument in order for this to work!
        (maybe sloppy? not sure how else to manage it other than established kwarg conventions)

        :param function:    some instance method that outputs a figure to the "outpath" argument.
        :param kwargs:      dictionary of any other keyword arguments for that function.
        :return:
        """

        if kwargs is None:
            kwargs = {}

        kwargs['outpath'] = self.absolute_filepath      # add outpath kwarg
        try:
            function(**kwargs)
        except TypeError:
            traceback.print_exc()
            raise TypeError("function {0} did not accept arguments {1}".format(function, kwargs))


    def get_tex(self):
        texdata = ["\\begin{figure}[H]",
                   "\\centering",

                   # removes extension from relative filepath
                   "\\includegraphics[width={width}]{{{rpath}}}".format(
                       width=self.width,
                       rpath="".join(self.relative_filepath.replace("\\", "/").split(".")[:-1])),

                   "\\caption{{{caption}}}".format(caption=self.caption),
                   "\\label{{{label}}}".format(label=self.label),
                   "\\end{figure}",
                   "\n"]
        return texdata
