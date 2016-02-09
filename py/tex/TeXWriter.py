__author__ = 'Jwely'

import os
from TeXFigureGenerator import TeXFigureGenerator


class TeXWriter:

    def __init__(self, main_path, texfile_path):

        self.main_path = os.path.abspath(main_path)
        self.texfile_path = os.path.abspath(texfile_path)

        self.content = []           # a list of content strings

        if not os.path.exists(self.main_path):
            raise Warning("main_path '{0}' does not exist".format(self.main_path))


    def add_figure(self, figure_path, caption, width,
                   create_from_function=None, create_kwargs=None):
        """
        Adds a figure to the TexFile, if create params are given, the figure is
        generated from scratch.
        :param figure_path:     filepath to the image
        :param caption:         caption to use on the figure
        :param width:           width of the figure as a string, such as '5in'

        :param create_from_function: input method to use to generate the figure
        :param create_kwargs: kwargs for the create_from_function
        """

        fig = TeXFigureGenerator(self.main_path, figure_path, caption, width)
        self.content += fig.get_tex()

        # if creation parameters were given, create the figure from scratch
        if create_from_function is not None:
            fig.create_from(create_from_function, create_kwargs)


    def add_text(self, text):
        """
        Adds a string of text to the TeX document, must be TeX formatted. Not really intended
        to be used extensively, but could be needed.
        """
        self.content += [text]


    def chapter(self, name):
        self.content += ["\\chapter{{{0}}}".format(name)]


    def write(self, verbose=False, reset_content=False, include_labels=True):
        """ Writes a `.tex` file with the current text and figures """

        with open(self.texfile_path, 'w+') as f:
            for line in self.content:
                if not include_labels:
                    if r"\label" not in line:
                        if verbose:
                            print(line)
                        f.write(line + "\n")
                else:
                    if verbose:
                        print(line)
                    f.write(line + "\n")
        print("Wrote tex file at {0}".format(self.texfile_path))

        # resets the content to ensure figure is only writen once, not subsequently
        if reset_content:
            self.content = []
