# -*- coding: utf-8 -*-
"""This module contains classes for calculating functional molecular information statistics.

"""
from collections import Counter, defaultdict
from multiprocessing import Pool
from operator import itemgetter
from string import Template
from ast import literal_eval as make_tuple
import bisect
import copy
import pkgutil
import itertools
import sys
import numpy as np
from pandas.io.formats.format import return_docstring

import nsb_entropy as nb
# import tsfm.nsb_entropy as nb
import random
import time
import math as mt
import exact as exact
# import tsfm.exact as exact
import glob
import re
import statsmodels.stats.multitest as smm
import pandas as pd


class DistanceCalculator:
    """A `DistanceCalculator` object contains methods for calculating several pairwise distance metrics between function logos.
    
    Currently, a `DistanceCalculator` object can calculate pairwise distance using the square-root of the Jensen-Shannon
    divergence and will print the resulting distance matrix to stdout.

    Args:
        distance (str): Indicates which distance metric to use for pairwise calculations.

    Attributes:
        distanceMetric (str): Indicates the distance metric to be used in 
            pairwise calculations.
        featureSet (:obj:`set` of :obj:`str`): A :obj:`set` of the structural 
            features contained in the function logos being compared (e.g. 1A, 173AU).
        functionSet (:obj:`set` of :obj:`str`): A :obj:`set` of the functional 
            classes contained in the function logos being compared.

    Example::
        
        x = tsfm.MolecularInformation.DistanceCalculator('jsd')
        x.get_distance(function_logos)

    """

    def __init__(self, distance):
        """The initialization of a `DistanceCalculator` object requires a :str: indicating the distance metric to be used.

        """
        self.distanceMetric = distance
        self.featureSet = set()
        self.functionSet = set()

    def get_distance(self, ResultsDict):
        """
        Prints a pairwise distance matrix using the distance metric indicated during instantiation to file. 

        Args:
            ResultsDict (:obj:`dict` of :obj:`str` mapping to :class:`FunctionLogoResults`):
                The values of the :obj:`dict` are compared using the selected pairwise
                distance metric.
        
        
        Note:
            Creates a :obj:`dict` of :obj:`str`: :class:`pandas.DataFrame` from 
            :obj:`ResultsDict`. The index of the dataframes are the union 
            of the structural features contained in :obj:`ResultsDict`, 
            and columns labels are the union of the functional classes contained in 
            :obj:`ResultsDict` including  a column containing 
            the functional information of the feature measured in bits. 
            Rows contain the Gorodkin fractional heights of each functional 
            class of each feature along with the functional information of the 
            feature measured in bits. The fractional heights of 
            each row is normalized to account for filtering of data and rounding 
            errors. The :obj:`dict` of :obj:`str`: :obj:`pandas.DataFrame` is 
            passed to the distance method set when the :class:`DistanceCalculator` 
            was instantiated. Below is an example of the :class:`pandas.DataFrame` 
            created\:
            
            +--------+-------+-------+-------+-------+-------+-------+--------+
            |        |   A   |   C   |   D   |   E   |   F   |   E   |  bits  |
            +========+=======+=======+=======+=======+=======+=======+========+
            |   1A   | 0.500 | 0.250 | 0.125 | 0.000 | 0.000 | 0.125 | 2.453  |
            +--------+-------+-------+-------+-------+-------+-------+--------+
            |   1U   | 0.000 | 0.250 | 0.125 | 0.500 | 0.125 | 0.000 | 2.453  |
            +--------+-------+-------+-------+-------+-------+-------+--------+

        """
        for result in ResultsDict:
            for coord in ResultsDict[result].basepairs:
                if (coord in ResultsDict[result].info):
                    for pairtype in ResultsDict[result].info[coord]:
                        self.featureSet.add("{}{}".format("".join(str(i) for i in coord), pairtype))
                        for function in ResultsDict[result].height[coord][pairtype]:
                            self.functionSet.add(function)

            for coord in range(ResultsDict[result].pos):
                if (coord in ResultsDict[result].info):
                    for base in ResultsDict[result].info[coord]:
                        self.featureSet.add("{}{}".format(coord, base))
                        for function in ResultsDict[result].height[coord][base]:
                            self.functionSet.add(function)

            # add inverse info features
            for coord in ResultsDict[result].basepairs:
                if (coord in ResultsDict[result].inverseInfo):
                    for pairtype in ResultsDict[result].inverseInfo[coord]:
                        self.featureSet.add("i{}{}".format("".join(str(i) for i in coord), pairtype))

            for coord in range(ResultsDict[result].pos):
                if (coord in ResultsDict[result].inverseInfo):
                    for base in ResultsDict[result].inverseInfo[coord]:
                        self.featureSet.add("i{}{}".format(coord, base))

        # remove features that contain gaps
        self.featureSet = {feature for feature in self.featureSet if not "-" in feature}

        # prepare pandas dataframes for each result object
        functionDict = {}
        pandasDict = {}
        for function in self.functionSet:
            functionDict[function] = np.zeros(len(self.featureSet), )
        functionDict["bits"] = np.zeros(len(self.featureSet), )

        for result in ResultsDict:
            pandasDict[result] = pd.DataFrame(functionDict, index=self.featureSet)
            for coord in ResultsDict[result].basepairs:
                if (coord in ResultsDict[result].info):
                    for pairtype in [pair for pair in ResultsDict[result].info[coord] if not "-" in pair]:
                        row = "{}{}".format("".join(str(i) for i in coord), pairtype)
                        pandasDict[result].loc[row, "bits"] = ResultsDict[result].info[coord][pairtype]
                        for function in ResultsDict[result].height[coord][pairtype]:
                            pandasDict[result].loc[row, function] = ResultsDict[result].height[coord][pairtype][
                                function]

            for coord in range(ResultsDict[result].pos):
                if (coord in ResultsDict[result].info):
                    for base in [nuc for nuc in ResultsDict[result].info[coord] if not nuc == "-"]:
                        row = "{}{}".format(coord, base)
                        pandasDict[result].loc[row, "bits"] = ResultsDict[result].info[coord][base]
                        for function in ResultsDict[result].height[coord][base]:
                            pandasDict[result].loc[row, function] = ResultsDict[result].height[coord][base][function]

            for coord in ResultsDict[result].basepairs:
                if (coord in ResultsDict[result].inverseInfo):
                    for pairtype in [pair for pair in ResultsDict[result].inverseInfo[coord] if not "-" in pair]:
                        row = "i{}{}".format("".join(str(i) for i in coord), pairtype)
                        pandasDict[result].loc[row, "bits"] = ResultsDict[result].inverseInfo[coord][pairtype]
                        for function in ResultsDict[result].inverseHeight[coord][pairtype]:
                            pandasDict[result].loc[row, function] = ResultsDict[result].inverseHeight[coord][pairtype][
                                function]

            for coord in range(ResultsDict[result].pos):
                if (coord in ResultsDict[result].inverseInfo):
                    for base in [nuc for nuc in ResultsDict[result].inverseInfo[coord] if not nuc == "-"]:
                        row = "i{}{}".format(coord, base)
                        pandasDict[result].loc[row, "bits"] = ResultsDict[result].inverseInfo[coord][base]
                        for function in ResultsDict[result].inverseHeight[coord][base]:
                            pandasDict[result].loc[row, function] = ResultsDict[result].inverseHeight[coord][base][
                                function]

        # normalize heights to equal one after possible removal of CIFs based on some criteria
        for frame in pandasDict:
            pandasDict[frame] = pandasDict[frame].round(3)
            pandasDict[frame].drop('bits', axis=1).div(pandasDict[frame].drop('bits', axis=1).sum(axis=1), axis=0)

        if (self.distanceMetric == "jsd"):
            self.rJSD(pandasDict)

    def rJSD(self, pandasDict):
        """
        Produces pairwise comparisons using rJSD metric

        This is method should not be directly called. Instead use the
        :meth:`get_distance`. All pairwise comparsions of OTUs are produced
        and :meth:`rJSD_distance` is called to do the calculations.

        Args:
            pandasDict (:obj:`dict` of `str` mapping to :class:`pandas.DataFrame`): 
                See :meth:`get_distance` for the format of the Data Frames.
        """
        pairwise_combinations = itertools.permutations(pandasDict.keys(), 2)
        jsdDistMatrix = pd.DataFrame(index=list(pandasDict.keys()), columns=list(pandasDict.keys()))
        jsdDistMatrix = jsdDistMatrix.fillna(0)
        for pair in pairwise_combinations:
            distance = 0
            for i, row in pandasDict[pair[0]].iterrows():
                if (row['bits'] == 0 and pandasDict[pair[1]].loc[i, 'bits'] == 0):
                    continue
                else:
                    distance += self.rJSD_distance(row.drop('bits').as_matrix(),
                                                   pandasDict[pair[1]].loc[i,].drop('bits').as_matrix(),
                                                   row['bits'], pandasDict[pair[1]].loc[i, 'bits'])

            jsdDistMatrix.loc[pair[0], pair[1]] = distance

        jsdDistMatrix = jsdDistMatrix.round(6)
        jsdDistMatrix.to_csv("jsdDistance.matrix", sep="\t")

    def entropy(self, dist):
        return np.sum(-dist[dist != 0] * np.log2(dist[dist != 0]))

    def rJSD_distance(self, dist1, dist2, pi1, pi2):
        r"""
        Weighted square root of the generalized Jensen-Shannon divergence defined by Lin 1991

        .. math::  
            
            D(X,Y) \equiv \sum_{f \in F} (I_f^X + I_f^Y) \sqrt{H[\pi_f^X p_f^X + \pi_f^Y p_f^Y] - (\pi_f^X H[p_f^X] + \pi_f^Y H[p_f^Y])}

        where :math:`\pi_f^X = \frac{I_f^X}{I_f^X + I_f^Y}` and :math:`\pi_f^Y = \frac{I_f^Y}{I_f^X + I_f^Y}`

        """
        step = self.entropy(pi1 * dist1 + pi2 * dist2) - (pi1 * self.entropy(dist1) + pi2 * self.entropy(dist2))
        return (pi1 + pi2) * mt.sqrt(step if step >= 0 else 0)


class FunctionLogoResults:
    """
    Stores results from information calculations and provides methods for text output and visualization.

    Args:
        name (:obj:`str`): Value is used as prefix for output files.
        basepairs (:obj:`list` of :obj:`tuples` of (:obj:`int`, :obj:`int`)):
            a list of basepair coordinates encoded as a :obj:`tuple` of two 
            :obj:`int`. 
            
            Note:
                This data structure is created as an attribute of 
                :class:`FunctionLogo` during instantiation and can be accessed
                with :attr:`FunctionLogo.basepairs` or created during
                instantiation of this class when ``from_file = True``
        pos (:obj:`int`): Stores length of the alignment.
            
            Note:
                See note for :attr:`basepairs`. Accessed using :attr:`FunctionLogo.pos`.
        sequences (:obj:`list` of :class:`Seq`): a list of :class:`Seq` objects 
            used for text output and visualization.

            Note:
                See note for :attr:`basepairs`. Accessed using :attr:`FunctionLogo.seq`
        pairs (:obj:`set` of :obj:`str`): unique basepair states found in the dataset.

            Note: 
                See note for :attr:`basepairs`.
        singles (:obj:`set` of :obj:`str`): unique states for single sites.
            
            Note:
                See note for :attr:`basepairs`.
        
        info (:obj:`dict` of :obj:`int` or :obj:`tuple` mapping to :obj:`dict` of :obj:`str` mapping to :obj:`float`):
            mapping of structural features to information content. Add this data structure using :meth:`add_information`.

            Note:
                This data structure is output of 
                :meth:`FunctionLogo.calculate_entropy_NSB()` or 
                :meth:`FunctionLogo.calculate_entropy_MM()`.
        height (:obj:`dict` of :obj:`int` or :obj:`tuple` mapping to :obj:`dict` of :obj:`str` mapping to :obj:`dict` of :obj:`str` mapping to :obj:`float`):
            mapping of structural features and functional class to class height. Add this data structure using :meth:`add_information`.

            Note:
                This data structure is output of 
                :meth:`FunctionLogo.calculate_entropy_NSB()` or 
                :meth:`FunctionLogo.calculate_entropy_MM()`.
        inverseInfo (:obj:`dict` of :obj:`int` or :obj:`tuple` mapping to :obj:`dict` of :obj:`str` mapping to :obj:`float`):
            mapping of structural features to information content for anti-determinants. Add this data structure using :meth:`add_information`.

            Note:
                This data structure is output of 
                :meth:`FunctionLogo.calculate_entropy_inverse_NSB()` or 
                :meth:`FunctionLogo.calculate_entropy_inverse_MM()`.
        inverseHeight (:obj:`dict` of :obj:`int` or :obj:`tuple` mapping to :obj:`dict` of :obj:`str` mapping to :obj:`dict` of :obj:`str` mapping to :obj:`float`):
            mapping of structural features and functional class to class height for anti-determinants. Add this data structure using :meth:`add_information`.

            Note:
                This data structure is output of 
                :meth:`FunctionLogo.calculate_entropy_inverse_NSB()` or 
                :meth:`FunctionLogo.calculate_entropy_inverse_MM()`.
        p (:obj:`dict` of :obj:`str` mapping to :obj:`dict`): mapping of structural features and class height to p-values.
            
            Note:
                This data structure is created using :meth:`add_stats()`
        inverse_p (:obj:`dict` of :obj:`str` mapping to :obj:`dict`): mapping of structural features and class height to p-values for anti-determinants
            
            Note:
                This data structure is created using :meth:`add_stats()`
        from_file (:obj:`bool`): create :class:`FunctionLogoResults` 
            object from file written with 
            :meth:`FunctionLogResults.text_output`
                
    """

    def __init__(self, name, basepairs=None, pos=0, sequences=None, pairs=None, singles=None, info=None,
                 height=None, inverseInfo=None, inverseHeight=None, p=None,
                 inverse_p=None, from_file=False):
        self.pos = pos
        self.correction = ""
        if (not info):
            self.info = defaultdict(lambda: defaultdict(float))
            self.height = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        else:
            self.info = info
            self.height = height

        if (not inverseInfo):
            self.inverseInfo = defaultdict(lambda: defaultdict(float))
            self.inverseHeight = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        else:
            self.inverseInfo = inverseInfo
            self.inverseHeight = inverseHeight

        if (not p):
            self.p = {'P': defaultdict(lambda: defaultdict(float)),
                      'p': defaultdict(lambda: defaultdict(lambda: defaultdict(float))),
                      'P_corrected': defaultdict(lambda: defaultdict(float)),
                      'p_corrected': defaultdict(lambda: defaultdict(lambda: defaultdict(float)))}
        else:
            self.p = p

        if (not inverse_p):
            self.inverse_p = {'P': defaultdict(lambda: defaultdict(float)),
                              'p': defaultdict(lambda: defaultdict(lambda: defaultdict(float))),
                              'P_corrected': defaultdict(lambda: defaultdict(float)),
                              'p_corrected': defaultdict(lambda: defaultdict(lambda: defaultdict(float)))}
        else:
            self.inverse_p = inverse_p

        if (not basepairs):
            self.basepairs = []
        else:
            self.basepairs = basepairs
        if (not sequences):
            self.sequences = []
        else:
            self.sequences = sequences

        if (not pairs):
            self.pairs = set()
        else:
            self.pairs = pairs

        if (not singles):
            self.singles = set()
        else:
            self.singles = singles

        if (from_file):
            self.name = name.split("/")[-1]
            self.from_file(name)
        else:
            self.name = name

    def from_file(self, file_name):
        """
        Read previously calculated results from file.

        Populates :class:`FunctionLogoResults` from previously calculated 
        results written to a file using :meth:`text_output`. 

        Args:
            file_name(:obj:`str`): File path of previously caclulated results
        """
        pvalue = False
        file_handle = open(file_name, "r")
        for line in file_handle:
            if (line.startswith("#")):
                if ("p-value" in line):
                    pvalue = True
            else:
                line = line.strip()
                spline = line.split("\t")
                if (spline[0] == "bp:"):
                    if (not make_tuple(spline[1]) in self.basepairs):
                        self.basepairs.append(make_tuple(spline[1]))
                    self.pairs.add(spline[2])
                    self.info[make_tuple(spline[1])][spline[2]] = float(spline[4])
                    if (pvalue):
                        self.p['P'][make_tuple(spline[1])][spline[2]] = float(spline[5])
                        self.p['P_corrected'][make_tuple(spline[1])][spline[2]] = float(spline[6])
                    for function in spline[7].split():
                        function_split = function.split(":")
                        self.height[make_tuple(spline[1])][spline[2]][function_split[0]] = float(function_split[1])
                        if (pvalue):
                            self.p['p'][make_tuple(spline[1])][spline[2]][function_split[0]] = float(function_split[2])
                            self.p['p_corrected'][make_tuple(spline[1])][spline[2]][function_split[0]] = float(
                                function_split[3])
                elif (spline[0] == "ss:"):
                    if (self.pos < int(spline[1])):
                        self.pos = int(spline[1])
                    self.singles.add(spline[2])
                    self.info[int(spline[1])][spline[2]] = float(spline[4])
                    if (pvalue):
                        self.p['P'][int(spline[1])][spline[2]] = float(spline[5])
                        self.p['P_corrected'][int(spline[1])][spline[2]] = float(spline[6])
                    for function in spline[7].split():
                        function_split = function.split(":")
                        self.height[int(spline[1])][spline[2]][function_split[0]] = float(function_split[1])
                        if (pvalue):
                            self.p['p'][int(spline[1])][spline[2]][function_split[0]] = float(function_split[2])
                            self.p['p_corrected'][int(spline[1])][spline[2]][function_split[0]] = float(
                                function_split[3])
                elif (spline[0] == "ibp:"):
                    if (not make_tuple(spline[1]) in self.basepairs):
                        self.basepairs.append(make_tuple(spline[1]))
                    self.pairs.add(spline[2])
                    self.inverseInfo[make_tuple(spline[1])][spline[2]] = float(spline[4])
                    if (pvalue):
                        self.inverse_p['P'][make_tuple(spline[1])][spline[2]] = float(spline[5])
                        self.inverse_p['P_corrected'][make_tuple(spline[1])][spline[2]] = float(spline[6])
                    for function in spline[7].split():
                        function_split = function.split(":")
                        self.inverseHeight[make_tuple(spline[1])][spline[2]][function_split[0]] = float(
                            function_split[1])
                        if (pvalue):
                            self.inverse_p['p'][make_tuple(spline[1])][spline[2]][function_split[0]] = float(
                                function_split[2])
                            self.inverse_p['p_corrected'][make_tuple(spline[1])][spline[2]][function_split[0]] = float(
                                function_split[3])
                elif (spline[0] == "iss:"):
                    if (self.pos < int(spline[1])):
                        self.pos = int(spline[1])
                    self.singles.add(spline[2])
                    self.inverseInfo[int(spline[1])][spline[2]] = float(spline[4])
                    if (pvalue):
                        self.inverse_p['P'][int(spline[1])][spline[2]] = float(spline[5])
                        self.inverse_p['P_corrected'][int(spline[1])][spline[2]] = float(spline[6])
                    for function in spline[7].split():
                        function_split = function.split(":")
                        self.inverseHeight[int(spline[1])][spline[2]][function_split[0]] = float(function_split[1])
                        if (pvalue):
                            self.inverse_p['p'][int(spline[1])][spline[2]][function_split[0]] = float(function_split[2])
                            self.inverse_p['p_corrected'][int(spline[1])][spline[2]][function_split[0]] = float(
                                function_split[3])
        self.pos += 1  # fix off by one
        file_handle.close()

    def add_information(self, info, height, inverse=False):
        """
        Add data structures containing results from information calculations

        This method is used to add results from 
        :meth:`FunctionLogo.calculate_entropy_NSB()`, 
        :meth:`FunctionLogo.calculate_entropy_MM()`,
        :meth:`FunctionLogo.calculate_entropy_inverse_NSB()` or 
        :meth:`FunctionLogo.calculate_entropy_inverse_MM()`. If reading previous
        results from a file this method is unnecessary because these data structures
        are populated from values in the file.

        Args:
            info (:obj:`dict`): mapping of structural features to information 
                content. This data structure is output of 
                :meth:`FunctionLogo.calculate_entropy_NSB()` or 
                :meth:`FunctionLogo.calculate_entropy_MM()`.
            height (:obj:`dict`): mapping of structural features and functional class to class height.
                This data structure is output of 
                :meth:`FunctionLogo.calculate_entropy_NSB()` or 
                :meth:`FunctionLogo.calculate_entropy_MM()`.
            inverse (:obj:`bool`): Defines if the data structures are for
                anti-determinates.
        """
        if (inverse):
            self.inverseInfo = info
            self.inverseHeight = height
        else:
            self.info = info
            self.height = height

    def add_stats(self, distribution, correction, inverse=False):
        """
        Perform statisical testing and multiple test correction

        Calculates p-values and multiple testing corrected p-values for
        structural features and functional class heights. Requires an
        instance of :class:`FunctionLogoDist` and calls the 
        :meth:`FunctionLogoDist.stat_test`. Methods for multiple test
        correction are provided by :class:`statsmodels.stats.multitest`.

        Args:
            distribution (:class:`FunctionLogoDist`): discrete probability 
                distributions of information content of structural 
                features and functional class height.
            correction (:obj:`str`): Multiple test correction method.
            inverse (:obj:`bool`): Produce statistical tests for
                anti-determinates.
        """
        self.correction = correction
        if (inverse):
            self.inverse_p = distribution.stat_test(self.inverseInfo, self.inverseHeight,
                                                    correction)
        else:
            self.p = distribution.stat_test(self.info, self.height, correction)

    def get(self, position, state):
        ret_counter = Counter()
        if (len(position) == 1):
            for x in self.sequences:
                if (x.seq[position[0]] == state[0]):
                    ret_counter[x.function] += 1
        if (len(position) == 2):
            for x in self.sequences:
                if (x.seq[position[0]] == state[0] and x.seq[position[1]] == state[1]):
                    ret_counter[x.function] += 1

        return ret_counter

    def text_output(self):
        """
        Write results to file named\: :attr:`name`\_results.txt
        """
        # build output heading
        file_handle = open("{}_results.txt".format(self.name.split("/")[-1]), "w")
        heading_dict = {}
        if (self.p):
            heading_dict['P'] = "\tp-value\t{}".format(self.correction)
            heading_dict['p'] = "\tclass:height:p-value:{}".format(self.correction)
        else:
            heading_dict['P'] = ""
            heading_dict['p'] = "\tclass:height"

        print("#bp\tcoord\tstate\tN\tinfo{P}{p}".format(**heading_dict), file=file_handle)
        for coord in sorted(self.basepairs, key=itemgetter(0)):
            if (coord in self.info):
                for pairtype in sorted(self.info[coord]):
                    output_string = "bp:\t{}".format(coord)
                    output_string += "\t{}\t{}\t{:05.3f}\t".format(pairtype, sum(self.get(coord, pairtype).values()),
                                                                   self.info[coord][pairtype])
                    if (self.p):
                        output_string += "{:08.6f}".format(self.p['P'][coord][pairtype])
                        output_string += "\t{:08.6f}".format(self.p['P_corrected'][coord][pairtype])

                    output_string += "\t"
                    for aainfo in sorted(self.height[coord][pairtype].items(), key=itemgetter(1), reverse=True):
                        output_string += "{}:{:05.3f}".format(aainfo[0], aainfo[1])
                        if (self.p):
                            output_string += ":{:08.6f}".format(self.p['p'][coord][pairtype][aainfo[0].upper()])
                            output_string += ":{:08.6f}".format(
                                self.p['p_corrected'][coord][pairtype][aainfo[0].upper()])
                        output_string += " "

                    print(output_string, file=file_handle)

        if (self.inverseInfo):
            print("#ibp\tcoord\tstate\tN\tinfo{P}{p}".format(**heading_dict), file=file_handle)
        for coord in sorted(self.basepairs, key=itemgetter(0)):
            if (coord in self.inverseInfo):
                for pairtype in sorted(self.inverseInfo[coord]):
                    output_string = "ibp:\t{}".format(coord)
                    output_string += "\t{}\t{}\t{:05.3f}\t".format(pairtype, sum(self.get(coord, pairtype).values()),
                                                                   self.inverseInfo[coord][pairtype])
                    if (self.p):
                        output_string += "{:08.6f}".format(self.inverse_p['P'][coord][pairtype])
                        output_string += "\t{:08.6f}".format(self.inverse_p['P_corrected'][coord][pairtype])

                    output_string += "\t"
                    for aainfo in sorted(self.inverseHeight[coord][pairtype].items(), key=itemgetter(1), reverse=True):
                        output_string += "{}:{:05.3f}".format(aainfo[0], aainfo[1])
                        if (self.p):
                            output_string += ":{:08.6f}".format(self.inverse_p['p'][coord][pairtype][aainfo[0].upper()])
                            output_string += ":{:08.6f}".format(
                                self.inverse_p['p_corrected'][coord][pairtype][aainfo[0].upper()])
                        output_string += " "

                    print(output_string, file=file_handle)

        print("#ss\tcoord\tstate\tN\tinfo{P}{p}".format(**heading_dict), file=file_handle)
        for coord in range(self.pos):
            if (coord in self.info):
                for base in sorted(self.info[coord]):
                    output_string = "ss:\t{}\t{}\t{}\t{:05.3f}".format(coord, base,
                                                                       sum(self.get([coord], base).values()),
                                                                       self.info[coord][base])
                    if (self.p):
                        output_string += "\t{:08.6f}".format(self.p['P'][coord][base])
                        output_string += "\t{:08.6f}".format(self.p['P_corrected'][coord][base])

                    output_string += "\t"
                    for aainfo in sorted(self.height[coord][base].items(), key=itemgetter(1), reverse=True):
                        output_string += "{}:{:05.3f}".format(aainfo[0], aainfo[1])
                        if (self.p):
                            output_string += ":{:08.6f}".format(self.p['p'][coord][base][aainfo[0].upper()])
                            output_string += ":{:08.6f}".format(self.p['p_corrected'][coord][base][aainfo[0].upper()])
                        output_string += " "

                    print(output_string, file=file_handle)

        if (self.inverseInfo):
            print("#iss\tcoord\tstate\tN\tinfo{P}{p}".format(**heading_dict), file=file_handle)
        for coord in range(self.pos):
            if (coord in self.inverseInfo):
                for base in sorted(self.inverseInfo[coord]):
                    output_string = "iss:\t{}\t{}\t{}\t{:05.3f}".format(coord, base,
                                                                        sum(self.get([coord], base).values()),
                                                                        self.inverseInfo[coord][base])
                    if (self.p):
                        output_string += "\t{:08.6f}".format(self.inverse_p['P'][coord][base])
                        output_string += "\t{:08.6f}".format(self.inverse_p['P_corrected'][coord][base])

                    output_string += "\t"
                    for aainfo in sorted(self.inverseHeight[coord][base].items(), key=itemgetter(1), reverse=True):
                        output_string += "{}:{:05.3f}".format(aainfo[0], aainfo[1])
                        if (self.p):
                            output_string += ":{:08.6f}".format(self.inverse_p['p'][coord][base][aainfo[0].upper()])
                            output_string += ":{:08.6f}".format(
                                self.inverse_p['p_corrected'][coord][base][aainfo[0].upper()])
                        output_string += " "

                    print(output_string, file=file_handle)
        file_handle.close()

    def logo_output_idLogo(self):
        coord_length = 0  # used to determine eps height
        coord_length_addition = 0
        logo_outputDict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

        for coord in range(self.pos):
            for base in sorted(self.singles):

                if (base in self.info[coord]):
                    if base not in self.height[coord]:
                        self.height[coord][base] = {}
                    for aainfo in sorted(self.height[coord][base].items(), key=itemgetter(1), reverse=True):
                        logo_outputDict[base][coord][aainfo[0]] = self.info[coord][base] * aainfo[1]
                else:
                    logo_outputDict[base][coord] = {}

        # output logos
        for base in logo_outputDict:
            logodata = ""
            for coord in sorted(logo_outputDict[base].keys()):
                if (len(str(coord)) > coord_length):
                    coord_length = len(str(coord))
                logodata += "numbering {{({}) makenumber}} if\ngsave\n".format(coord)
                for aainfo in sorted(logo_outputDict[base][coord].items(), key=itemgetter(1)):
                    if (aainfo[1] < 0.0001 or mt.isnan(aainfo[1])):
                        continue
                    logodata += "{:07.5f} ({}) numchar\n".format(aainfo[1], aainfo[0].upper())
                logodata += "grestore\nshift\n"
            # output logodata to template
            template_byte = pkgutil.get_data('tsfm', 'eps/Template.eps')
            logo_template = template_byte.decode('utf-8')
            with open("{}_{}.eps".format(base, self.name.split("/")[-1]), "w") as logo_output:
                src = Template(logo_template)
                if (len(base) == 2):
                    logodata_dict = {'logo_data': logodata, 'low': min(logo_outputDict[base].keys()),
                                     'high': max(logo_outputDict[base].keys()),
                                     'length': 21 * len(logo_outputDict[base].keys()),
                                     'height': 735 - (5 * (coord_length + coord_length_addition))}
                else:
                    logodata_dict = {'logo_data': logodata, 'low': min(logo_outputDict[base].keys()),
                                     'high': max(logo_outputDict[base].keys()),
                                     'length': 15.68 * len(logo_outputDict[base].keys()),
                                     'height': 735 - (5 * (coord_length + coord_length_addition))}
                logo_output.write(src.substitute(logodata_dict))

    def logo_output(self, inverse=False):
        """
        Produce function logo postscript files
        """
        coord_length = 0  # used to determine eps height
        coord_length_addition = 0
        # print(" hiiiiiiiiiiiii", len(self.basepairs))
        logo_outputDict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

        # logo output dict construction
        for coord in sorted(self.basepairs, key=itemgetter(0)):

            for pairtype in sorted(self.pairs):
                if (coord in self.info and pairtype in self.info[coord]):  # we added this.
                    for aainfo in sorted(self.height[coord][pairtype].items(), key=itemgetter(1), reverse=True):
                        logo_outputDict[pairtype][coord][aainfo[0]] = self.info[coord][pairtype] * aainfo[1]
                else:
                    logo_outputDict[pairtype][coord] = {}

        for coord in range(self.pos):  # 0,1,2..,72
            for base in sorted(self.singles):  # ACGU..
                if (base in self.info[coord]):
                    # print(self.height[coord][base].items(), " items.")
                    for aainfo in sorted(self.height[coord][base].items(), key=itemgetter(1), reverse=True):
                        logo_outputDict[base][coord][aainfo[0]] = self.info[coord][base] * aainfo[1]
                else:
                    logo_outputDict[base][coord] = {}

        # output logos
        for base in logo_outputDict:
            logodata = ""
            for coord in sorted(logo_outputDict[base].keys()):
                if (len(str(coord)) > coord_length):
                    coord_length = len(str(coord))
                logodata += "numbering {{({}) makenumber}} if\ngsave\n".format(coord)
                for aainfo in sorted(logo_outputDict[base][coord].items(), key=itemgetter(1)):
                    if (aainfo[1] < 0.0001 or mt.isnan(aainfo[1])):
                        continue
                    logodata += "{:07.5f} ({}) numchar\n".format(aainfo[1], aainfo[0].upper())
                logodata += "grestore\nshift\n"
            # output logodata to template
            template_byte = pkgutil.get_data('tsfm', 'eps/Template.eps')
            logo_template = template_byte.decode('utf-8')
            with open("{}_{}.eps".format(base, self.name.split("/")[-1]), "w") as logo_output:
                src = Template(logo_template)
                if (len(base) == 2):
                    logodata_dict = {'logo_data': logodata, 'low': min(logo_outputDict[base].keys()),
                                     'high': max(logo_outputDict[base].keys()),
                                     'length': 21 * len(logo_outputDict[base].keys()),
                                     'height': 735 - (5 * (coord_length + coord_length_addition))}
                    # print(735 - (5 * (coord_length + coord_length_addition)), "heigh base")
                else:
                    logodata_dict = {'logo_data': logodata, 'low': min(logo_outputDict[base].keys()),
                                     'high': max(logo_outputDict[base].keys()),
                                     'length': 15.68 * len(logo_outputDict[base].keys()),
                                     'height': 735 - (5 * (coord_length + coord_length_addition))}
                    # print(735 - (5 * (coord_length + coord_length_addition)), "heigh", base)

                    # print(logo_outputDict[base].keys(),base)
                    # print(logodata_dict)
                    # print( 735 - (5 * (coord_length + coord_length_addition)),"kkkkkkkkkkkkkk")
                logo_output.write(src.substitute(logodata_dict))

        if (inverse):

            inverse_logo_outputDict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
            # inverse logo output dict construction
            for coord in sorted(self.basepairs, key=itemgetter(0)):
                for pairtype in sorted(self.pairs):
                    if (pairtype in self.inverseInfo[coord]):
                        for aainfo in sorted(self.inverseHeight[coord][pairtype].items(), key=itemgetter(1),
                                             reverse=True):
                            inverse_logo_outputDict[pairtype][coord][aainfo[0]] = self.inverseInfo[coord][pairtype] * \
                                                                                  aainfo[1]
                    else:
                        inverse_logo_outputDict[pairtype][coord] = {}

            for coord in range(self.pos):
                for base in sorted(self.singles):
                    if (base in self.inverseInfo[coord]):
                        for aainfo in sorted(self.inverseHeight[coord][base].items(), key=itemgetter(1), reverse=True):
                            inverse_logo_outputDict[base][coord][aainfo[0]] = self.inverseInfo[coord][base] * aainfo[1]
                    else:
                        inverse_logo_outputDict[base][coord] = {}

            for base in inverse_logo_outputDict:
                logodata = ""
                for coord in sorted(inverse_logo_outputDict[base].keys()):
                    if (len(str(coord)) > coord_length):
                        coord_length = len(str(coord))
                    logodata += "numbering {{({}) makenumber}} if\ngsave\n".format(coord)
                    for aainfo in sorted(inverse_logo_outputDict[base][coord].items(), key=itemgetter(1)):
                        if (aainfo[1] < 0.0001 or mt.isnan(aainfo[1])):
                            continue
                        logodata += "{:07.5f} ({}) numchar\n".format(aainfo[1], aainfo[0].upper())
                    logodata += "grestore\nshift\n"
                # output logodata to template
                template_byte = pkgutil.get_data('tsfm', 'eps/Template.eps')
                logo_template = template_byte.decode('utf-8')
                with open("inverse_{}_{}.eps".format(base, self.name.split("/")[-1]), "w") as logo_output:
                    src = Template(logo_template)
                    if (len(base) == 2):
                        logodata_dict = {'logo_data': logodata, 'low': min(inverse_logo_outputDict[base].keys()),
                                         'high': max(inverse_logo_outputDict[base].keys()),
                                         'length': 21 * len(inverse_logo_outputDict[base].keys()),
                                         'height': 735 - (5 * (coord_length + coord_length_addition))}
                    else:
                        logodata_dict = {'logo_data': logodata, 'low': min(inverse_logo_outputDict[base].keys()),
                                         'high': max(inverse_logo_outputDict[base].keys()),
                                         'length': 15.68 * len(inverse_logo_outputDict[base].keys()),
                                         'height': 735 - (5 * (coord_length + coord_length_addition))}
                    logo_output.write(src.substitute(logodata_dict))


class FunctionLogoDist:
    """
    Discrete probability distributions of information values.
    Probabilty distributions are created using a permutation label shuffling
    strategy. Permuted data is created using :meth:`FunctionLogo.permute` and
    distribution are inferred from the permuted data and 
    :class:`FunctionLogoDist` objects created using 
    :meth:`FunctionLogo.permInfo`.

    Attributes:
        bpinfodist (:obj:`dict` of :obj:`float` mapping to :obj:`int`):
            Discrete probability distribution of basepair feature information
        bpheightdist (:obj:`dict` of :obj:`float` mapping to :obj:`int`):
            Discrete probability distribution of functional class 
            information of basepair features
        singleinfodist (:obj:`dict` of :obj:`float` mapping to :obj:`int`):
            Discrete probability distribution of single base feature information
        singleheightdist (:obj:`dict` of :obj:`float` mapping to :obj:`int`):
            Discrete probability distribution of functional class 
            information of single base features

    """

    def __init__(self):

        self.bpinfodist = defaultdict(int)
        self.bpheightdist = defaultdict(int)

        self.singleinfodist = defaultdict(int)
        self.singleheightdist = defaultdict(int)

    def weighted_dist(self, bpdata, singledata):
        for x in bpdata[0]:
            self.bpinfodist[x] += 1
        for x in bpdata[1]:
            self.bpheightdist[x] += 1

        for x in singledata[0]:
            self.singleinfodist[x] += 1
        for x in singledata[1]:
            self.singleheightdist[x] += 1

        self.bpinfo_sorted_keys = sorted(self.bpinfodist.keys())
        self.bpheight_sorted_keys = sorted(self.bpheightdist.keys())
        self.ssinfo_sorted_keys = sorted(self.singleinfodist.keys())
        self.ssheight_sorted_keys = sorted(self.singleheightdist.keys())

    def stat_test(self, info, height, correction):
        """
        Performs statistical tests and multiple test correction.

        Calculates a p-value using a right tail probability test on the
        instance's discrete probability distributions. Methods for multiple test
        correction are provided by :class:`statsmodels.stats.multitest`. This
        method is usually invoked using :meth:`FunctionLogoResults.add_stats`.

        Args:
            info (:obj:`dict` of :obj:`int` or :obj:`tuple` mapping to :obj:`dict` of :obj:`str` mapping to :obj:`float`):
                mapping of structural features to information content.
            height (:obj:`dict` of :obj:`int` or :obj:`tuple` mapping to :obj:`dict` of :obj:`str` mapping to :obj:`dict` of :obj:`str` mapping to :obj:`float`):
                mapping of structural features and functional class to class height.
            correction (:obj:`str`): Method for multiple test correction. Any 
                method available in :class:`statsmodels.stats.multitest` is a
                valid option
        """
        P = defaultdict(lambda: defaultdict(float))
        P_corrected = defaultdict(lambda: defaultdict(float))
        p = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        p_corrected = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        bp_coords = []
        ss_coords = []
        for coord in info:
            for pairtype in info[coord]:
                if ("," in str(coord)):
                    bp_coords.append(coord)
                    P[coord][pairtype] = self.rtp(self.bpinfodist, info[coord][pairtype], self.bpinfo_sorted_keys)
                    for aa in height[coord][pairtype]:
                        p[coord][pairtype][aa] = self.rtp(self.bpheightdist,
                                                          info[coord][pairtype] * height[coord][pairtype][aa],
                                                          self.bpheight_sorted_keys)
                else:
                    ss_coords.append(coord)
                    P[coord][pairtype] = self.rtp(self.singleinfodist, info[coord][pairtype], self.ssinfo_sorted_keys)
                    for aa in height[coord][pairtype]:
                        p[coord][pairtype][aa] = self.rtp(self.singleheightdist,
                                                          info[coord][pairtype] * height[coord][pairtype][aa],
                                                          self.ssheight_sorted_keys)
        test_bp = []
        test_ss = []
        bp_coords.sort()
        ss_coords.sort()
        for coord in bp_coords:
            for pairtype in sorted(P[coord]):
                test_bp.append(P[coord][pairtype])

        for coord in ss_coords:
            for pairtype in sorted(P[coord]):
                test_ss.append(P[coord][pairtype])

        test_bp_results = smm.multipletests(test_bp, method=correction)[1].tolist()
        test_ss_results = smm.multipletests(test_ss, method=correction)[1].tolist()

        for coord in bp_coords:
            for pairtype in sorted(P[coord]):
                P_corrected[coord][pairtype] = test_bp_results.pop(0)

        for coord in ss_coords:
            for pairtype in sorted(P[coord]):
                P_corrected[coord][pairtype] = test_ss_results.pop(0)

        test_bp = []
        test_ss = []
        for coord in bp_coords:
            for pairtype in sorted(p[coord]):
                for aa in sorted(p[coord][pairtype]):
                    test_bp.append(p[coord][pairtype][aa])

        for coord in ss_coords:
            for pairtype in sorted(p[coord]):
                for aa in sorted(p[coord][pairtype]):
                    test_ss.append(p[coord][pairtype][aa])

        test_bp_results = smm.multipletests(test_bp, method=correction)[1].tolist()
        test_ss_results = smm.multipletests(test_ss, method=correction)[1].tolist()

        for coord in bp_coords:
            for pairtype in sorted(p[coord]):
                for aa in sorted(p[coord][pairtype]):
                    p_corrected[coord][pairtype][aa] = test_bp_results.pop(0)

        for coord in ss_coords:
            for pairtype in sorted(p[coord]):
                for aa in sorted(p[coord][pairtype]):
                    p_corrected[coord][pairtype][aa] = test_ss_results.pop(0)

        return {'P': P, 'p': p, "P_corrected": P_corrected, "p_corrected": p_corrected}

    def rtp(self, data, point, keys_sorted):
        if (point > 0):
            part = 0
            total = sum(data.values())
            i = bisect.bisect_left(keys_sorted, point)
            if (point <= keys_sorted[-1]):
                for y in keys_sorted[i:]:
                    part += data[y]
                return part / total
            else:
                return 0.0
        else:
            return 1.0


class Seq:
    """
    Providing a data structure constisting of a molecular sequence labeled with a functional class.

    Args:
        function (:obj:`str`): Functional annotation of the sequence.
        seq (:obj:`str`): Molecular sequence data.
    """

    def __init__(self, function, seq):
        self.function = function
        self.seq = seq

    def __len__(self):
        return len(self.seq)


class FunctionLogo:
    """
    Parses structural and sequence infomation and provides methods for Function Logo calculations

    This class provided data structures and methods for calculating 
    functional information of basepair a single base features. Additionally,
    methods for producing permuted data sets with function class labels
    shuffled.
    
    Args:
        struct_file (:obj:`str`): File name containing secondary structure
            notation in cove, infernal, or text format.
        kind (:obj:`str`): secondary structure notation format.

    """

    def __init__(self, struct_file, kind=None, exact_init=None, inverse_init=None):
        if (exact_init):
            self.exact = exact_init
        else:
            self.exact = []

        if (inverse_init):
            self.inverse_exact = inverse_init
        else:
            self.inverse_exact = []

        if (kind):
            if (kind == "s"):
                self.basepairs = []
            else:
                self.parse_struct(struct_file, kind)
        else:
            self.basepairs = struct_file
        self.pos = 0
        self.sequences = []
        self.pairs = set()
        self.singles = set()
        self.functions = Counter()

    def parse_sequences(self, file_prefix):
        """
        Parse sequence alignment data in clustal format

        Sequence alignment files are required to be in clustal format with
        each functional class having its own file. Alignment files must
        conform to the naming standard ``fileprefix_functionalclass.aln``.

        Args:
            file_prefix (:obj:`str`): Prefix used to identify a group of alignment files.

        """
        for fn in glob.glob("{}_?.aln".format(file_prefix)):
            match = re.search("_([A-Z])\.aln", fn)
            aa_class = match.group(1)
            with open(fn, "r") as ALN:
                good = False
                begin_seq = False
                interleaved = False
                seq = {}
                for line in ALN:
                    match = re.search("^(\S+)\s+(\S+)", line)
                    if (re.search("^CLUSTAL", line)):
                        good = True
                        continue
                    elif (re.search("^[\s\*\.\:]+$", line) and not interleaved and begin_seq):
                        interleaved = True
                    elif (re.search("^[\s\*\.\:]+$", line) and interleaved and begin_seq):
                        continue
                    elif (match and not interleaved):
                        begin_seq = True
                        if (not good):
                            sys.exit("File {} appears not to be a clustal file".format(fn))
                        seq[match.group(1)] = match.group(2)
                    elif (match and interleaved):
                        seq[match.group(1)] += match.group(2)
            for sequence in seq.values():
                # print(aa_class)
                self.add_sequence(aa_class, sequence.upper().replace("T", "U"))

        print("{} alignments parsed".format(len(self.functions.keys())), file=sys.stderr)

    def parse_struct(self, struct_file, kind):
        """
        Parse secondary structure file for basepair locations.

        Args:
        struct_file (:obj:`str`): File containing structural annotation
        kind (:obj:`str`): Structural annotation format
        """
        print("Parsing base-pair coordinates", file=sys.stderr)
        basepairs = []
        ss = ""
        pairs = defaultdict(list)
        tarm = 0
        stack = []
        if (kind == "infernal"):
            for line in struct_file:
                line = line.strip()
                ss += line.split()[2]
            struct_file.seek(0)

            state = "start"
            for count, i in enumerate(ss):
                if (i == "("):
                    if (state == "start"):
                        state = "A"
                elif (i == "<"):
                    stack.append(count)
                    if (state == "A"):
                        state == "D"
                    elif (state == "cD"):
                        state = "C"
                    elif (state == "cC"):
                        state = "T"
                elif (i == ">"):
                    if (state == "D"):
                        state = "cD"
                    elif (state == "C"):
                        state = "cC"
                    elif (state == "T"):
                        state = "cT"

                    arm = state.replace("c", "")
                    pairs[arm].append([stack.pop(), count])
                elif (i == ")"):
                    pairs['A'].append([stack.pop(), count])

            for arm in pairs:
                for pair in pairs[arm]:
                    basepairs.append((pair[0], pair[1]))

        if (kind == "cove"):
            for line in struct_file:
                line = line.strip()
                ss += line.split()[1]
            struct_file.seek(0)

            state = "start"
            for count, i in enumerate(ss):
                if (i == ">" and (state == "start" or state == "AD")):
                    if (state == "start"):
                        state = "AD"
                    stack.append(count)

                elif (i == "<" and (state == "AD" or state == "D")):
                    if (state == "AD"):
                        state = "D"
                    pairs[state].append([stack.pop(), count])

                elif (i == ">" and (state == "D" or state == "C")):
                    if (state == "D"):
                        state = "C"
                    stack.append(count)

                elif (i == "<" and (state == "C" or state == "cC")):
                    if (state == "C"):
                        state = "cC"
                    pairs["C"].append([stack.pop(), count])

                elif (i == ">" and (state == "cC" or state == "T")):
                    if (state == "cC"):
                        state = "T"
                    stack.append(count)
                    tarm += 1

                elif (i == "<" and (state == "T" and tarm > 0)):
                    pairs[state].append([stack.pop(), count])
                    tarm -= 1

                elif (i == "<" and (state == "T" or state == "A") and tarm == 0):
                    state = "A"
                    pairs[state].append([stack.pop(), count])

            for arm in pairs:
                for pair in pairs[arm]:
                    basepairs.append((pair[0], pair[1]))

        if (kind == "text"):
            for line in struct_file:
                coords = "".join(line.split(":")[1])
                coords = coords.split(",")
                for coord1, coord2 in zip(coords[0::2], coords[1::2]):
                    basepairs.append((int(coord1), int(coord2)))

        self.basepairs = basepairs

    def approx_expect(self, H, k, N):
        return H - ((k - 1) / ((mt.log(4)) * N))

    def exact_run(self, n, p, numclasses):
        ######################################
        j = exact.calc_exact(n, p, numclasses)
        print("{:2} {:07.5f}".format(n, j[1]), file=sys.stderr)
        return j

    def permuted(self, items, pieces=2):
        random.seed()
        sublists = [[] for i in range(pieces)]
        for x in items:
            sublists[random.randint(0, pieces - 1)].append(x)
        permutedList = []
        for i in range(pieces):
            time.sleep(0.01)
            random.seed()
            random.shuffle(sublists[i])
            permutedList.extend(sublists[i])
        return permutedList

    def permutations(self, numPerm, aa_classes):
        indices = []
        permStructList = []
        for p in range(numPerm):
            indices.append(self.permuted(aa_classes))
        for index in indices:
            permStruct = FunctionLogo(self.basepairs, exact_init=self.exact, inverse_init=self.inverse_exact)
            for i, seqs in enumerate(self.sequences):
                permStruct.add_sequence(index[i], seqs.seq)
            permStructList.append(permStruct)
        return permStructList

    def permute(self, permute_num, proc):
        """
        Creates permuted datasets by shuffling functional annotation labels of sequences.

        Args:
            permute_num (:obj:`int`): Number of permutations to perform
            proc (:obj:`int`): Number of concurrent processes to run
        """
        with Pool(processes=proc) as pool:
            perm_jobs = []
            for x in range(proc):
                if (x == 0):
                    perm_jobs.append((permute_num // proc + permute_num % proc, self.get_functions()))
                else:
                    perm_jobs.append((permute_num // proc, self.get_functions()))

            perm_results = pool.starmap(self.permutations, perm_jobs)
            self.permutationList = []
            for x in perm_results:
                self.permutationList += x

    def permInfo(self, method, proc, inverse=False):
        """
        Calculate functional information statistics of permuted datasets.

        Args:
            method (:obj:`str`): Entropy estimation method. Either NSB or Miller-Maddow.
            proc (:obj:`int`): Number of concurrent processes to run.

        Return:
        perm_dist (:class:`FunctionLogoDist`): Discrete distribution of 
            functional information estimated from permuted datasets.
        """
        bp_info = []
        bp_height = []
        single_info = []
        single_height = []
        with Pool(processes=proc) as pool:
            if (len(self.permutationList) < proc):
                chunk = 1
            else:
                chunk = len(self.permutationList) // proc

            if (not inverse):
                if (method == "NSB"):
                    perm_info_results = pool.map(self.perm_info_calc_NSB, self.permutationList, chunk)
                else:
                    perm_info_results = pool.map(self.perm_info_calc_MM, self.permutationList, chunk)
            else:
                if (method == "NSB"):
                    perm_info_results = pool.map(self.perm_info_calc_inverse_NSB, self.permutationList, chunk)
                else:
                    perm_info_results = pool.map(self.perm_info_calc_inverse_MM, self.permutationList, chunk)

        for perm in perm_info_results:
            bp_info.extend(perm[0])
            single_info.extend(perm[1])
            bp_height.extend(perm[2])
            single_height.extend(perm[3])

        perm_dist = FunctionLogoDist()
        perm_dist.weighted_dist((bp_info, bp_height), (single_info, single_height))
        return perm_dist

    def perm_info_calc_MM(self, x):
        total_info_bp = []
        height_info_bp = []
        total_info_ss = []
        height_info_ss = []
        info, height_dict = x.calculate_entropy_MM()
        for coord in sorted(self.basepairs, key=itemgetter(0)):
            if (coord in info):
                for pairtype in sorted(info[coord]):
                    total_info_bp.append(info[coord][pairtype])
                    for aainfo in sorted(height_dict[coord][pairtype].items(), key=itemgetter(1), reverse=True):
                        height_info_bp.append(aainfo[1] * info[coord][pairtype])

        for coord in range(self.pos):
            if (coord in info):
                for base in sorted(info[coord]):
                    total_info_ss.append(info[coord][base])
                    for aainfo in sorted(height_dict[coord][base].items(), key=itemgetter(1), reverse=True):
                        height_info_ss.append(aainfo[1] * info[coord][base])

        return (total_info_bp, total_info_ss, height_info_bp, height_info_ss)

    def perm_info_calc_inverse_MM(self, x):
        total_info_bp = []
        height_info_bp = []
        total_info_ss = []
        height_info_ss = []
        info, height_dict = x.calculate_entropy_inverse_MM()

        for coord in sorted(self.basepairs, key=itemgetter(0)):
            if (coord in info):
                for pairtype in sorted(info[coord]):
                    total_info_bp.append(info[coord][pairtype])
                    for aainfo in sorted(height_dict[coord][pairtype].items(), key=itemgetter(1), reverse=True):
                        height_info_bp.append(aainfo[1] * info[coord][pairtype])

        for coord in range(self.pos):
            if (coord in info):
                for base in sorted(info[coord]):
                    total_info_ss.append(info[coord][base])
                    for aainfo in sorted(height_dict[coord][base].items(), key=itemgetter(1), reverse=True):
                        height_info_ss.append(aainfo[1] * info[coord][base])

        return (total_info_bp, total_info_ss, height_info_bp, height_info_ss)

    def perm_info_calc_inverse_NSB(self, x):
        total_info_bp = []
        height_info_bp = []
        total_info_ss = []
        height_info_ss = []
        info, height_dict = x.calculate_entropy_inverse_NSB()

        for coord in sorted(self.basepairs, key=itemgetter(0)):
            if (coord in info):
                for pairtype in sorted(info[coord]):
                    total_info_bp.append(info[coord][pairtype])
                    for aainfo in sorted(height_dict[coord][pairtype].items(), key=itemgetter(1), reverse=True):
                        height_info_bp.append(aainfo[1] * info[coord][pairtype])

        for coord in range(self.pos):
            if (coord in info):
                for base in sorted(info[coord]):
                    total_info_ss.append(info[coord][base])
                    for aainfo in sorted(height_dict[coord][base].items(), key=itemgetter(1), reverse=True):
                        height_info_ss.append(aainfo[1] * info[coord][base])

        return (total_info_bp, total_info_ss, height_info_bp, height_info_ss)

    def perm_info_calc_NSB(self, x):
        total_info_bp = []
        height_info_bp = []
        total_info_ss = []
        height_info_ss = []
        info, height_dict = x.calculate_entropy_NSB()

        for coord in sorted(self.basepairs, key=itemgetter(0)):
            if (coord in info):
                for pairtype in sorted(info[coord]):
                    total_info_bp.append(info[coord][pairtype])
                    for aainfo in sorted(height_dict[coord][pairtype].items(), key=itemgetter(1), reverse=True):
                        height_info_bp.append(aainfo[1] * info[coord][pairtype])

        for coord in range(self.pos):
            if (coord in info):
                for base in sorted(info[coord]):
                    total_info_ss.append(info[coord][base])
                    for aainfo in sorted(height_dict[coord][base].items(), key=itemgetter(1), reverse=True):
                        height_info_ss.append(aainfo[1] * info[coord][base])

        return (total_info_bp, total_info_ss, height_info_bp, height_info_ss)

    def calculate_exact(self, n, proc, inverse=False):
        #
        # number of CPUs in the system.
        """
        Exact method of small sample size correction.

        Calculate the exact method of sample size correction for up to N samples.
        Computational intensive portion of the calculation is implemented as a C
        extension. This method is fully described in Schneider et al 1986. 
        This calculation is polynomial in sample size. It becomes prohibitively 
        expensive to calculate beyond a sample size of 16. The correction 
        factor of each sample size will be calculated in parallel up to 
        :obj:`proc` at a time.

        Args:
            n (:obj:`int`): Calculate correction up to this sample size.
            proc (:obj:`int`): Number of concurrent processes to run
            inverse (:obj:`bool`): If true calculate sample size correction
                for anti-determinates.
        """
        exact_list = []
        exact_results = []
        if (inverse):
            inverse_functions = Counter()
            for aa_class in self.functions:
                inverse_functions[aa_class] = sum(self.functions.values()) / self.functions[aa_class]

            p = [x / sum(list(inverse_functions.values())) for x in inverse_functions.values()]
            for i in range(1, n + 1):
                exact_list.append((i, p, len(self.functions.values())))

            with Pool(processes=proc) as pool:
                exact_results = pool.starmap(self.exact_run, exact_list)

            for x in exact_results:
                self.inverse_exact.append(x[1])
        else:
            p = [x / sum(list(self.functions.values())) for x in
                 self.functions.values()]  # PORTION of each function key. for example: x:0.34
            for i in range(1, n + 1):
                exact_list.append((i, p, len(self.functions.values())))

            with Pool(processes=proc) as pool:
                exact_results = pool.starmap(self.exact_run, exact_list)
                # print(exact_list,"exact_list................")

            for x in exact_results:
                self.exact.append(x[1])
                # exact is an array of length n = 10

    def calculate_entropy_MM(self):
        """
        Calculate functional information using Miller-Maddow estimator.
        """
        info = defaultdict(lambda: defaultdict(float))
        height_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

        functions_array = np.array(list(self.functions.values()))
        bg_entropy = -np.sum(
            (functions_array[functions_array != 0] / functions_array[functions_array != 0].sum()) * np.log2(
                functions_array[functions_array != 0] / functions_array[functions_array != 0].sum()))
        for pairs in self.basepairs:
            for state in self.pairs:
                state_counts = self.get(pairs, state)
                if (sum(state_counts.values()) == 0):
                    continue

                nsb_array = np.array(list(state_counts.values()) + [0] * (len(self.functions) - len(state_counts)))
                fg_entropy = -np.sum((nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()) * np.log2(
                    nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()))
                if (sum(state_counts.values()) <= len(self.exact)):
                    expected_bg_entropy = self.exact[sum(state_counts.values()) - 1]
                else:
                    expected_bg_entropy = self.approx_expect(bg_entropy, len(self.functions),
                                                             sum(state_counts.values()))

                if (expected_bg_entropy - fg_entropy < 0):
                    info[pairs][state] = 0
                else:
                    info[pairs][state] = expected_bg_entropy - fg_entropy

                height_class = {}
                for aa_class in state_counts:
                    height_class[aa_class] = (state_counts[aa_class] / sum(state_counts.values())) / (
                            self.functions[aa_class] / len(self))
                for aa_class in height_class:
                    height_dict[pairs][state][aa_class] = height_class[aa_class] / sum(height_class.values())

        for singles in range(self.pos):
            for state in self.singles:
                state_counts = self.get([singles], state)
                if (sum(state_counts.values()) == 0):
                    continue

                nsb_array = np.array(list(state_counts.values()) + [0] * (len(self.functions) - len(state_counts)))
                fg_entropy = -np.sum((nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()) * np.log2(
                    nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()))
                if (sum(state_counts.values()) <= len(self.exact)):
                    expected_bg_entropy = self.exact[sum(state_counts.values()) - 1]
                else:
                    expected_bg_entropy = self.approx_expect(bg_entropy, len(self.functions),
                                                             sum(state_counts.values()))

                if (expected_bg_entropy - fg_entropy < 0):
                    info[singles][state] = 0
                else:
                    info[singles][state] = expected_bg_entropy - fg_entropy

                height_class = {}
                for aa_class in state_counts:
                    height_class[aa_class] = (state_counts[aa_class] / sum(state_counts.values())) / (
                            self.functions[aa_class] / len(self))
                for aa_class in height_class:
                    height_dict[singles][state][aa_class] = height_class[aa_class] / sum(height_class.values())

        return (info, height_dict)

    def calculate_entropy_inverse_MM(self):
        """
        Calculate functional information for anit-determinates using Miller-Maddow estimator.
        """
        info_inverse = defaultdict(lambda: defaultdict(float))
        height_dict_inverse = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        inverse_functions = Counter()
        for aa_class in self.functions:
            inverse_functions[aa_class] = sum(self.functions.values()) / self.functions[aa_class]

        np_inverse_functions = np.array(list(inverse_functions.values()))
        bg_entropy = -np.sum((np_inverse_functions[np_inverse_functions != 0] / np_inverse_functions[
            np_inverse_functions != 0].sum()) * np.log2(
            np_inverse_functions[np_inverse_functions != 0] / np_inverse_functions[np_inverse_functions != 0].sum()))
        for pairs in self.basepairs:
            for state in self.pairs:
                state_counts = self.get(pairs, state)
                if (sum(state_counts.values()) == 0):
                    continue
                if (not len(state_counts) == len(self.functions)):
                    for function in self.functions:
                        state_counts[function] += 1

                inverse_state_counts = Counter()
                for aa_class in state_counts:
                    inverse_state_counts[aa_class] = sum(state_counts.values()) / state_counts[aa_class]

                nsb_array = np.array(list(inverse_state_counts.values()))
                fg_entropy = -np.sum((nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()) * np.log2(
                    nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()))
                if (sum(state_counts.values()) <= len(self.inverse_exact)):
                    expected_bg_entropy = self.inverse_exact[sum(state_counts.values()) - 1]
                else:
                    expected_bg_entropy = self.approx_expect(bg_entropy, len(self.functions),
                                                             sum(state_counts.values()))

                if (expected_bg_entropy - fg_entropy < 0):
                    info_inverse[pairs][state] = 0
                else:
                    info_inverse[pairs][state] = expected_bg_entropy - fg_entropy

                height_class = {}
                for aa_class in inverse_state_counts:
                    height_class[aa_class] = (inverse_state_counts[aa_class] / sum(inverse_state_counts.values())) / (
                            inverse_functions[aa_class] / sum(inverse_functions.values()))
                for aa_class in height_class:
                    height_dict_inverse[pairs][state][aa_class] = height_class[aa_class] / sum(height_class.values())

        for singles in range(self.pos):
            for state in self.singles:
                state_counts = self.get([singles], state)
                if (sum(state_counts.values()) == 0):
                    continue
                if (not len(state_counts) == len(self.functions)):
                    for function in self.functions:
                        state_counts[function] += 1

                inverse_state_counts = Counter()
                for aa_class in state_counts:
                    inverse_state_counts[aa_class] = sum(state_counts.values()) / state_counts[aa_class]

                nsb_array = np.array(list(inverse_state_counts.values()))
                fg_entropy = -np.sum((nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()) * np.log2(
                    nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()))
                if (sum(state_counts.values()) <= len(self.inverse_exact)):
                    expected_bg_entropy = self.inverse_exact[sum(state_counts.values()) - 1]
                else:
                    expected_bg_entropy = self.approx_expect(bg_entropy, len(self.functions),
                                                             sum(state_counts.values()))

                if (expected_bg_entropy - fg_entropy < 0):
                    info_inverse[singles][state] = 0
                else:
                    info_inverse[singles][state] = expected_bg_entropy - fg_entropy

                height_class = {}
                for aa_class in inverse_state_counts:
                    height_class[aa_class] = (inverse_state_counts[aa_class] / sum(inverse_state_counts.values())) / (
                            inverse_functions[aa_class] / sum(inverse_functions.values()))
                for aa_class in height_class:
                    height_dict_inverse[singles][state][aa_class] = height_class[aa_class] / sum(height_class.values())

        return (info_inverse, height_dict_inverse)

    def calculate_entropy_inverse_NSB(self):
        """
        Calculate functional information for anit-determinates using NSB estimator.
        """
        # for singles in range(self.pos):
        #     for state in self.singles:
        #         state_counts = self.get([singles], state)
        #         for p in state_counts:
        #             print(state_counts[p],"state_counts[p]")
        #             print(state_counts[p].values(),"state_counts[p]")
        #             if state_counts[p] == 0:
        #                 print("hmmmmmmm")

        info_inverse = defaultdict(lambda: defaultdict(float))
        height_dict_inverse = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        inverse_functions = Counter()
        for aa_class in self.functions:
            inverse_functions[aa_class] = sum(self.functions.values()) / self.functions[aa_class]

        np_inverse_functions = np.array(list(inverse_functions.values()))
        bg_entropy = -np.sum((np_inverse_functions[np_inverse_functions != 0] / np_inverse_functions[
            np_inverse_functions != 0].sum()) * np.log2(
            np_inverse_functions[np_inverse_functions != 0] / np_inverse_functions[np_inverse_functions != 0].sum()))
        for pairs in self.basepairs:
            for state in self.pairs:
                state_counts = self.get(pairs, state)
                if (sum(state_counts.values()) == 0):
                    continue
                if (not len(state_counts) == len(self.functions)):
                    for function in self.functions:
                        state_counts[function] += 1

                inverse_state_counts = Counter()
                for aa_class in state_counts:
                    inverse_state_counts[aa_class] = sum(state_counts.values()) / state_counts[aa_class]

                nsb_array = np.array(list(inverse_state_counts.values()))
                if (sum(state_counts.values()) <= len(self.inverse_exact)):
                    expected_bg_entropy = self.inverse_exact[sum(state_counts.values()) - 1]
                    fg_entropy = -np.sum((nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()) * np.log2(
                        nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()))
                else:
                    expected_bg_entropy = bg_entropy
                    fg_entropy = nb.S(nb.make_nxkx(nsb_array, nsb_array.size), nsb_array.sum(), nsb_array.size)

                if (expected_bg_entropy - fg_entropy < 0):
                    info_inverse[pairs][state] = 0
                else:
                    info_inverse[pairs][state] = expected_bg_entropy - fg_entropy

                height_class = {}
                for aa_class in inverse_state_counts:
                    height_class[aa_class] = (inverse_state_counts[aa_class] / sum(inverse_state_counts.values())) / (
                            inverse_functions[aa_class] / sum(inverse_functions.values()))
                for aa_class in height_class:
                    height_dict_inverse[pairs][state][aa_class] = height_class[aa_class] / sum(height_class.values())

        for singles in range(self.pos):
            for state in self.singles:
                state_counts = self.get([singles], state)
                if (sum(state_counts.values()) == 0):
                    continue
                if (not len(state_counts) == len(self.functions)):
                    for function in self.functions:
                        state_counts[function] += 1

                inverse_state_counts = Counter()
                for aa_class in state_counts:
                    inverse_state_counts[aa_class] = sum(state_counts.values()) / state_counts[aa_class]

                nsb_array = np.array(list(inverse_state_counts.values()))
                if (sum(state_counts.values()) <= len(self.inverse_exact)):
                    expected_bg_entropy = self.inverse_exact[sum(state_counts.values()) - 1]
                    fg_entropy = -np.sum((nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()) * np.log2(
                        nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()))
                else:
                    expected_bg_entropy = bg_entropy
                    fg_entropy = nb.S(nb.make_nxkx(nsb_array, nsb_array.size), nsb_array.sum(), nsb_array.size)

                if (expected_bg_entropy - fg_entropy < 0):
                    info_inverse[singles][state] = 0
                else:
                    info_inverse[singles][state] = expected_bg_entropy - fg_entropy

                height_class = {}
                for aa_class in inverse_state_counts:
                    height_class[aa_class] = (inverse_state_counts[aa_class] / sum(inverse_state_counts.values())) / (
                            inverse_functions[aa_class] / sum(inverse_functions.values()))
                for aa_class in height_class:
                    height_dict_inverse[singles][state][aa_class] = height_class[aa_class] / sum(height_class.values())

        return (info_inverse, height_dict_inverse)

    def calculate_entropy_NSB(self):
        """
        Calculate functional information using NSB estimator.
        """
        info = defaultdict(lambda: defaultdict(float))
        height_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

        functions_array = np.array(list(self.functions.values()))

        # print(functions_array,"function_array")
        # print(self.functions,"self.function")
        # print(self.functions.values(),"self.functions.values")
        bg_entropy = -np.sum(
            (functions_array[functions_array != 0] / functions_array[functions_array != 0].sum()) * np.log2(
                functions_array[functions_array != 0] / functions_array[functions_array != 0].sum()))
        for pairs in self.basepairs:
            for state in self.pairs:
                state_counts = self.get(pairs, state)
                if (sum(state_counts.values()) == 0):
                    continue
                nsb_array = np.array(list(state_counts.values()) + [0] * (len(self.functions) - len(state_counts)))
                if (sum(state_counts.values()) <= len(self.exact)):
                    expected_bg_entropy = self.exact[sum(state_counts.values()) - 1]
                    fg_entropy = -np.sum((nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()) * np.log2(
                        nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()))
                else:
                    expected_bg_entropy = bg_entropy
                    fg_entropy = nb.S(nb.make_nxkx(nsb_array, nsb_array.size), nsb_array.sum(), nsb_array.size)

                if (expected_bg_entropy - fg_entropy < 0):
                    info[pairs][state] = 0
                else:
                    info[pairs][state] = expected_bg_entropy - fg_entropy

                height_class = {}
                for aa_class in state_counts:
                    height_class[aa_class] = (state_counts[aa_class] / sum(state_counts.values())) / (
                            self.functions[aa_class] / len(self))
                for aa_class in height_class:
                    height_dict[pairs][state][aa_class] = height_class[aa_class] / sum(height_class.values())

        for singles in range(self.pos):
            for state in self.singles:
                state_counts = self.get([singles], state)
                # print(state_counts,"statecounts")
                if (sum(state_counts.values()) == 0):
                    continue
                nsb_array = np.array(list(state_counts.values()) + [0] * (len(self.functions) - len(state_counts)))
                # print(self.exact,"self.exact",len(self.exact))
                if (sum(state_counts.values()) <= len(self.exact)):
                    expected_bg_entropy = self.exact[sum(state_counts.values()) - 1]
                    fg_entropy = -np.sum((nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()) * np.log2(
                        nsb_array[nsb_array != 0] / nsb_array[nsb_array != 0].sum()))
                else:
                    expected_bg_entropy = bg_entropy
                    fg_entropy = nb.S(nb.make_nxkx(nsb_array, nsb_array.size), nsb_array.sum(), nsb_array.size)

                if (expected_bg_entropy - fg_entropy < 0):
                    info[singles][state] = 0
                else:
                    info[singles][state] = expected_bg_entropy - fg_entropy

                height_class = {}
                for aa_class in state_counts:
                    height_class[aa_class] = (state_counts[aa_class] / sum(state_counts.values())) / (
                            self.functions[aa_class] / len(self))

                for aa_class in height_class:
                    height_dict[singles][state][aa_class] = height_class[aa_class] / sum(height_class.values())

        return (info, height_dict)

    def is_overlap(self, position):
        pass

    def add_sequence(self, function, seq):
        self.sequences.append(Seq(function, seq))
        # print(Seq(function, seq),"dddddddddddddddddddddd")
        self.functions[function] += 1
        self.pos = len(seq)  # pos is equal to 72
        self.singles.update(seq)  # this contains A U C G - N

        for x in self.basepairs:
            self.pairs.add(seq[x[0]] + seq[x[1]])

    def get(self, position, state):  # state is A U C G -
        ret_counter = Counter()

        if (len(position) == 1):
            for x in self.sequences:
                if (x.seq[position[0]] == state[0]):
                    ret_counter[x.function] += 1
        if (len(position) == 2):
            for x in self.sequences:
                if (x.seq[position[0]] == state[0] and x.seq[position[1]] == state[1]):
                    ret_counter[x.function] += 1
        return ret_counter

    def get_functions(self):
        function_list = []
        for key, val in self.functions.items():
            function_list.extend([key] * val)
        return function_list

    def __len__(self):
        return len(self.sequences)

    #  _______________________ ID logo Calculations ___________________________________________________

    def calculate_logoID_infos(self, info_1, info_2, pos, singles, pairs, basepairs):
        """
        Calculate information for id-logo (information difference of foreground and background).
        In id_info12, 2 is the background, 1 is foreground (ID_logo_info12 = info_1 - info_2)
        """

        id_info12 = defaultdict(lambda: defaultdict(float))
        id_info21 = defaultdict(lambda: defaultdict(float))
        for k in range(pos):

            logo_1 = info_1[k]
            logo_2 = info_2[k]

            for c in singles:
                id_info21[k][c] = logo_2[c] - logo_1[c]
                if id_info21[k][c] < 0:
                    id_info21[k][c] = 0

                id_info12[k][c] = logo_1[c] - logo_2[c]
                if id_info12[k][c] < 0:
                    id_info12[k][c] = 0

        for k in basepairs:

            logo_1 = info_1[k]
            logo_2 = info_2[k]
            for c in pairs:

                id_info21[k][c] = logo_2[c] - logo_1[c]
                if id_info21[k][c] < 0:
                    id_info21[k][c] = 0

                id_info12[k][c] = logo_1[c] - logo_2[c]
                if id_info12[k][c] < 0:
                    id_info12[k][c] = 0

        return id_info12, id_info21

    def calculate_logoID_heights(self, height_b, height_f, info, pos, singles, pairs, basepairs, type):
        """
        Calculate height of each symbol within a stack for id-logo.
        """

        id_height = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        # adding zero to all the functions that do not exist within a stack

        for single in range(pos):
            for state in singles:
                for p in height_f[single][state]:
                    if p in height_b[single][state] and height_b[single][state][p] != 0:
                        id_height[single][state][p] = (height_f[single][state][p] /
                                                       height_b[single][state][p])
                    else:
                        id_height[single][state][p] = height_f[single][state][p]
                summ = sum(id_height[single][state].values())
                for p in id_height[single][state]:
                    if summ != 0:
                        id_height[single][state][p] = id_height[single][state][p] * info[single][state] / summ

        for pair in basepairs:
            for state in pairs:
                for p in height_f[pair][state]:
                    if p in height_b[pair][state] and height_b[pair][state][p] != 0:
                        id_height[pair][state][p] = (height_f[pair][state][p] /
                                                     height_b[pair][state][p])
                    else:
                        id_height[pair][state][p] = height_f[pair][state][p]
                summ = sum(id_height[pair][state].values())
                for p in id_height[pair][state]:
                    if summ != 0:
                        id_height[pair][state][p] = id_height[pair][state][p] * info[pair][state] / summ

        for single in range(pos):
            for state in singles:
                mysum = sum(id_height[single][state].values())
                for p in id_height[single][state]:
                    if mysum != 0:
                        id_height[single][state][p] = id_height[single][state][p] / mysum

        for pair in basepairs:
            for state in pairs:
                mysum = sum(id_height[pair][state].values())
                for p in id_height[pair][state]:
                    if mysum != 0:
                        id_height[pair][state][p] = id_height[pair][state][p] / mysum

        # adding zero to all the functions that do not exist within a stack
        for single in basepairs:
            for state in pairs:
                for t in type:
                    if t not in id_height[single][state]:
                        id_height[single][state][t] = 0

        for single in range(pos):
            for state in singles:
                for t in type:
                    if t not in id_height[single][state]:
                        id_height[single][state][t] = 0

        return id_height

    # _______________________ KLD Difference logo Calculations ___________________________________________________

    def calculate_prob_dist(self, types, pairs, basepairs):
        """
        Calculate prior probability, p(y), for every function and posterior probability p(y|x) of each symbol within
        a stack for KLD-logo.
        This function will be called two times for the background and foreground
        """

        kld_post_dist = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        kld_prior_dist = defaultdict(float)
        functions_array = np.array(list(self.functions.values()))

        # adding one to the frequency of every function within a stack to avoid devide by zero!
        for p in types:
            if p not in self.functions:
                kld_prior_dist[p] = 1 / (functions_array[functions_array != 0].sum() + len(types))
            else:
                kld_prior_dist[p] = (self.functions[p] + 1) / (functions_array[functions_array != 0].sum() + len(types))

        for singles in range(self.pos):
            for state in self.singles:
                state_counts = self.get([singles], state)
                for t in types:
                    if t not in state_counts:
                        state_counts[t] = 1
                    else:
                        state_counts[t] += 1
                for p in state_counts:
                    kld_post_dist[singles][state][p] = state_counts[p] / sum(state_counts.values())

        for pair in basepairs:
            for state in pairs:
                state_counts = self.get(pair, state)
                for t in types:
                    if t not in state_counts:
                        state_counts[t] = 1
                    else:
                        state_counts[t] += 1

                for p in state_counts:
                    kld_post_dist[pair][state][p] = state_counts[p] / sum(state_counts.values())

        return kld_prior_dist, kld_post_dist

    def calculate_kld(self, back_prior, fore_prior, back_post, fore_post, pairs, basepairs):
        """
        Calculate height of each symbol within a stack for kld-logo.
        """

        # The height of the individual letters in a stack will be proportional to this ratio
        # ratio: a dictionary for keeping the ratios for each stack
        ratios = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

        kld_prior = 0

        # kld_dic: a dictionary for keeping the height of each stack in kld logo
        kld_dic = defaultdict(lambda: defaultdict(float))

        for p in back_prior:
            kld_prior += fore_prior[p] * np.log2(fore_prior[p] / back_prior[p])

        for singles in range(self.pos):
            for state in self.singles:
                for p in back_post[singles][state]:
                    kld_dic[singles][state] += fore_post[singles][state][p] * np.log2(
                        fore_post[singles][state][p] / back_post[singles][state][p])

                    ratios[singles][state][p] = (fore_post[singles][state][p] / fore_prior[p]) / (
                            back_post[singles][state][p] / back_prior[p])
                # the prior needs to be subtracted and the e(n) needs to be added when the sample size with a stack is
                # less than n (default 16)
                # kld_dic[singles][state] -= kld_prior

        for pair in basepairs:
            for state in pairs:
                for p in back_post[pair][state]:
                    kld_dic[pair][state] += fore_post[pair][state][p] * np.log2(
                        fore_post[pair][state][p] / back_post[pair][state][p])

                    ratios[pair][state][p] = (fore_post[pair][state][p] / fore_prior[p]) / (
                            back_post[pair][state][p] / back_prior[p])
                # kld_dic[pair][state] -= kld_prior

        # kld_heights: a dictionary for the height of each symbol within a stack for kld-logo.
        kld_heights = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

        for single in range(self.pos):
            for state in self.singles:
                for p in back_post[single][state]:
                    kld_heights[single][state][p] = kld_dic[single][state] * ratios[single][state][p] / sum(
                        ratios[single][state].values())

        for pair in basepairs:
            for state in pairs:
                for p in back_post[pair][state]:
                    kld_heights[pair][state][p] = kld_dic[pair][state] * ratios[pair][state][p] / sum(
                        ratios[pair][state].values())

        for single in range(self.pos):
            for state in self.singles:
                mysum = sum(kld_heights[single][state].values())
                for p in kld_heights[single][state]:
                    if mysum != 0:
                        kld_heights[single][state][p] = kld_heights[single][state][p] / mysum

        for pair in basepairs:
            for state in pairs:
                mysum = sum(kld_heights[pair][state].values())
                for p in kld_heights[pair][state]:
                    if mysum != 0:
                        kld_heights[pair][state][p] = kld_heights[pair][state][p] / mysum

        return kld_dic, kld_heights
