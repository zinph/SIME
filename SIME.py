#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      zinph
#
# Created:     01/11/2016
# Copyright:   (c) zinph 2016
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os
import re
import os.path
from molvs import standardize_smiles
from random import *
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from itertools import *
from collections import Counter, defaultdict
from rdkit.Chem.Scaffolds import MurckoScaffold
from operator import itemgetter

class SIME:

    def __init__(self):

##        self.info = {}  # to record number of compounds for each length
        self.total_numcompounds = 0
        self.load_sugars()
        self.load_extenders()
        self.library_size       = int(input('Desired Library Size (numbers only) :'))

##        self.SM_info = {}
        self.smile_file_name    = self.create_directory()+'/mcrl'
        self.max_repeat         = int(input('Maximum occurrence of the same structural motifs per scaffold (number only) :'))
        self.min_sugar          = int(input('Minimal number of sugars per scaffold (number only) :'))
        self.ext_question = input('Generate all possible stereocenters for extender structural motifs at joining carbons? y or n :')
        self.sugar_question = input('Generate all possible stereocenters for sugars at joining carbons? y or n :')


        self.info_manager = open(self.smile_file_name + '_info','a+')
        self.info_manager.write('Desired Library Size (numbers only) :'+ self.library_size + '\n')
        self.info_manager.write('Maximum occurrence of the same structural motifs per scaffold (number only) :'+ str(self.max_repeat) + '\n')
        self.info_manager.write('Minimal number of sugars per scaffold (number only) :'+ str(self.min_sugar) + '\n')
        self.info_manager.write('Generate all possible stereocenters for extender structural motifs at joining carbons? y or n :y'+ self.ext_question + '\n')
        self.info_manager.write('Generate all possible stereocenters for sugars at joining carbons? y or n :' + self.sugar_question + '\n')

    def create_directory(self):
        old_directory = os.getcwd()
        newfolder = input('Peferred Directory Name for Output Files: ')
        new_directory = os.path.join(old_directory, newfolder)
        while os.path.exists(new_directory):
            print('This folder exists or input is invalid. Try again.')
            newfolder = input('Peferred Directory Name for Output Files: ')
            new_directory = os.path.join(old_directory,newfolder)
        os.mkdir(new_directory)
        return new_directory


    def load_sugars(self):
        f = open('Data/sugars','r')
        original_sugars = f.read().splitlines()
        if self.sugar_question.lower() == 'no' or self.sugar_question.lower() == 'n':
            self.sugars = [r.replace('[*R*]','') for r in original_sugars]
        else:
            sugars = []
            for i in original_sugars:
                sugars.append(self.ENUMERATE_sugar_stereocenters(i))
            self.sugars = [r.replace('[*R*]','') for r in list(chain(*sugars))]
        self.make_full_sugar_list() # make self.full_list by adding hydroxyl to self.sugars
        self.info_manager.write('\n\nSugars\n'+'\n'.join(original_sugars)+'\n')

    def load_extenders(self):
        f = open('Data/selected_extenders.txt','r')
        original_extenders = f.read().splitlines()
        if self.ext_question.lower() == 'no' or self.ext_question.lower() == 'n':
             self.extenders = [r.replace('[*R*]','') for r in original_extenders]
        else:
            self.extenders = [r.replace('[*R*]','') for r in self.enumerate_SM_stereocenters(original_extenders)]
        self.info_manager.write('\n\nStructural Motifs\n'+'\n'.join(original_extenders)+'\n')
        self.info_manager.close()


    def make_full_sugar_list(self):
        '''
        self.full_list contains all the sugars and hydroxyl groups.
        '''
        hydroxyl = ['[C@H](O)','[C@@H](O)']
        self.full_list = self.sugars.copy()
        self.full_list += hydroxyl # contains all sugars and hydroxyl groups


    def ENUMERATE_sugar_stereocenters(self, smile):
        '''
        Take in sugar strings that start and end with [*R*], and return a list of sugars with two different stereoceters for the joining carbon.
        '''

        sugar_stereocenters = []
        # if the stereocenter of the joining carbon isn't defined
        if smile[5] is 'C':
#            template = smile[0:5] + '[C@H]' + smile[6:]
            template = smile.replace(smile[5], "[C@@H]", 1)
            sugar_stereocenters.append(template)
            template = smile.replace(smile[5], "[C@H]", 1)
            sugar_stereocenters.append(template)
        else:
            sugar_stereocenters.append(smile)
            if "@" in smile[:10]:
                if "@@" in smile[:10]:   # for clockwise
                    template = smile.replace('@@', '@', 1)
                else:
                    template = smile.replace('@', '@@', 1)
                sugar_stereocenters.append(template)
        return sugar_stereocenters


    def locate_SM_replace_points(self, smile):
        '''
        Take a string, and locate places for replacement. They are indicated by [1*], [2*], etc.... Return the string with all these joints replaced with [*]s.
        '''
        numbers = set(re.findall(r'\d+', smile))
        possible_joints = ['['+str(m) +'*]' for m in numbers]
        for each in possible_joints:
            smile = smile.replace(each, '[*]')
        return smile

    def generate_templates_withextenders(self, smile):
        '''
        Generate all possible templates. Takes in a smile string (structural core). This function only deals with extenders or structural motifs.
        Then, insert all possible extenders at those joint positions.
        '''
        smile_with_stars = [[r] for r in self.locate_SM_replace_points(smile).split('[*]')]  # take a string with [*]s and split into different fragments, convert each fragment into a list
        counter = 1
        # smile_with_stars = [ fragment1, fragment2, fragment3 ,...] all are split at joint positions
        # insert self.extenders in between all fragments (except for the first and last blocks).
        # so it will be something like [fragment1, [self.extenders], fragment2, [self.extenders], fragment3, [self.extenders] ,...]
        for i in range(len(smile_with_stars)-1):
            shuffle(self.extenders)
            smile_with_stars.insert(counter,self.extenders)
            counter+=2

        template = [x for x in smile_with_stars if x != ['']]
        self.make_compounds(template)
        return template


    def generate_templates_withExtendersNSugars(self,smile):
        '''

        '''
        smile_with_stars = self.string_splitter(self.locate_SM_replace_points(smile), '[*]')
        template = []
        # create a template holder that will have a list of extenders or sugars at the respective split location points and the rest of the core will remain the same.
        # The position of all these fragments (core, extenders, sugars) have to be in the correct order.
        # Each of all these fragments have to be in lists so you can perform product of them later.

        for each in range(len(smile_with_stars)):  # iterate for all the fragments initially split at extender/SM location points, sugar locations are embedded in some fragments within the list
            # works only for the already splitted list based on SM points '[*]'
            if '[*sugar*]' in smile_with_stars[each][0]: #locate sugar portions -- ['[*sugar*]'].  length of this is one and it should contain '[*sugar*]'
                sugar_fragments = self.string_splitter(smile_with_stars[each][0],'[*sugar*]')
                template+= sugar_fragments      # add sugar fragments
            else:
                template.append(smile_with_stars[each])

        template = [x for x in template if x != ['']]  # [stable_fragment1, [possible extender motifs], stable_fragment2, [possible sugar moieties], ...]
        SM_template = self.insert_SMs(template)  # a list of possible extenders inserted at SM locations
        SGR_order = self.generate_dummy_sugar_templates(SM_template,minimal_sugars=self.min_sugar)  # At least how many sugars do you want in the macrolide scaffold? because of this, more complications.

        for each in SGR_order:
            current_SYMBOLsugar_template = self.add_SYMBOLsugars_to_dummy_templates(each,SM_template) # add the lists of sugars and full_list at the dummy positions
            current_sugar_template = self.insert_sugars_to_dummies(current_SYMBOLsugar_template)
            self.make_compounds(current_sugar_template)



    def make_compounds(self,template):
        written = []
        file_counter = 1
        file_temp = self.smile_file_name + '_'+str(file_counter)+'.smiles' # attempts to split files because they get too large. Name of the first file will be "file_" + this variable
        file_handler = open(file_temp,'a+')
        for item in product(*template):
            if self.max_occurrence(list(item))[1] <= self.max_repeat: # If the count of most common SM is less than or equal to the number set up by the user

                if self.total_numcompounds <= self.library_size:
                    if len(written) < 1000000:
                        temp = ''.join([str(r) for r in item])
                        m = Chem.MolFromSmiles(temp)
                        self.total_numcompounds += 1
                        written.append(temp)
                    else:
                        file_handler.write('\n'.join(written))  # write smiles in written list
                        file_handler.close()
                        file_counter +=1
                        file_temp = self.smile_file_name + '_'+str(file_counter)+'.smiles'
                        file_handler = open(file_temp,'a+')
                        written = []
                else:
                    break
                    file_handler.close()
        file_handler.close()


    def RS_check(self,smile,ringsize):
        '''
        Take in a smile string and the desired ring size. For many of the macrolides, it will be 14.
        If the smile has the ring size of the desired number, return the same smile. If not, return None.
        '''
        m  = Chem.MolFromSmiles(smile)
        if m.GetAtomWithIdx(2).IsInRingSize(ringsize):
            return smile

    def max_occurrence(self, template_list):
        '''
        Take in a list of fragments, and return the most common fragment along with the number of occurrences.
        Returns a tuple of the most common fragment in SMILE format, and its occurrences.
        E.g.('CCCC',5)
        '''
#        most_common,num_most_common = Counter(template_list).most_common(1)[0] # SM, number of occurrences
        c = defaultdict(int)
        for item in template_list:
            c[item] += 1
        return max(c.items(), key=itemgetter(1))


    def enumerate_SM_stereocenters(self, ext):
        '''
        Enumerate possible stereostereocenters in the form of nested lists.
        For example, [[SM1_R, SM1_S],[SM2_R, SM2_S], ... , [plain SMs]]
        Append all other stereocenters of SMs in [R,S] and the rest of the plain SMs in one list at the end of the all_possible list.
        Then, combinations will be performed on this list to enumerate all possible templates.
        It returns a nested list templates with all possible stereocenters + one last list of plain SMs.

        handles SMs without any stereocenters like ketone
        '''
        SM_list = list(set(ext))
        all_stereos = []  # to store SMs with R and S stereocenters
        plain = []  # for SMs without any specified stereocenters
        for i in range(len(SM_list)):
            if "@" in SM_list[i]:
                if "@@" in SM_list[i]:   # for clockwise
                    new_stereo = SM_list[i].replace('@@','@',1)
                else:
                    new_stereo = SM_list[i].replace('@','@@',1)
                all_stereos.append([SM_list[i],new_stereo])   # Both stereocenters added as a list [R,S]
            else:
                plain.append(SM_list[i])  # for anticlockwise

        all_possible = all_stereos+[plain]
        all_possible = list(chain(*all_possible))

        return all_possible

    def replace_string_with_list(self,string, old, new):
        '''
        input  ---> string = '1[*]234[*]5[*]6', old = '[*]', new = ['a','b','c']
        output ---> [['1'], ['a', 'b', 'c'], ['234'], ['a', 'b', 'c'], ['5'], ['a', 'b', 'c'], ['6']]
        Takes in a string, split them at "old" positions at which "new" list is added.
        Returns a nested list of all fragments as shown in output.
        '''
        counter = 1
        frag = [[s] for s in string.split(old)]
        for i in range(len(frag)-1):
            frag.insert(counter,new)
            counter+=2
        return frag


    def string_splitter(self, string, split_character):
        '''
        Very similar to self.replace_string_with_list.
        input  ---> string = '1[*]234[*]5[*]6', old = '[*]'
        output ---> ['1', '[*]', '234', '[*]', '5', '[*]', '6']
        Takes in a string and split it at the positions of split character.
        Returns a list of all the split fragments, including the split characters. Each individual item in the list is a string.
        '''
        counter = 1
        frag = [[s] for s in string.split(split_character)]
        for i in range(len(frag)-1):
            frag.insert(counter,[split_character])
            counter+=2
        return frag

    def insert_SMs(self, template):
        '''
        Takes in a template resulted from string_splitter; it would look like
        --->
        template =[['1'], '[*]', ['2'],['[*sugar*]'],['3'], '[*]', ['[*sugar*]'],['4'], '[*]', ['5'], '[*]', ['6']]
        Replace ['[*]'] with SM list.
        '''
        # for SMs
        SM_indexes = [index for index, value in enumerate(template) if value == ['[*]']]
        for i in SM_indexes:
            template[i] = self.extenders
        # for sugars
        return template

    def generate_dummy_sugar_templates(self, template, minimal_sugars=1):
        '''
        n = the least number of sugars the users want in each macrolide.
        Default is one, i.e. there will be at least one sugar in each macrolide.
        Generate a list of all possible templates using dummys as 'SUGARS' (intended for only sugars) and 'FULL_LIST (intended for sugars + hydroxy).'

        '''
        num_sugars = template.count(['[*sugar*]'])
        list_with_atLeast_nSugars = n*['SUGARS']+(num_sugars-n)*['FULL_LIST']  # Make a new list with at least "n" "SUGARS" and the rest "FULL_LIST"
        sugar_lists_in_order = [] # make a new list to hold all possible sugar templates at each position
        # now create all different arrangements of sugars, full_list. The positions of these blocks matter.
        for i in permutations(list_with_atLeast_nSugars):
            if i not in sugar_lists_in_order:
                sugar_lists_in_order.append(i)
        return sugar_lists_in_order


    def add_SYMBOLsugars_to_dummy_templates(self, sugar_dummy_order,template_with_sugarinlist):
        '''
        each sugar_dummy_order looks like ('SUGARS', 'FULL_LIST', 'FULL_LIST', 'FULL_LIST')
        template_with_sugarinlist looks like = [['1'], [ext1,ext2,...], ['2'],['[*sugar*]'],['3'], [ext1,ext2,...], ['[*sugar*]'],['4'], [ext1,ext2,...], ['5'], [ext1,ext2,...], ['6']]
        '''
        copy = template_with_sugarinlist.copy()  # create a copy of the template_with_sugarinlist so it doesn't go and modify the original template
        sugar_indexes = [index for index, value in enumerate(template_with_sugarinlist) if value == ['[*sugar*]']]
        for i in range(len(sugar_indexes)):
            copy[sugar_indexes[i]] = sugar_dummy_order[i]  # match the dummies at the correct sugar positions
        return copy

    def insert_sugars_to_dummies(self, template):
        '''
        Takes in a template resulted from string_splitter; it would look like
        --->
        template =[['1'], '[*]', ['2'],['[*sugar*]'],['3'], '[*]', ['[*sugar*]'],['4'], '[*]', ['5'], '[*]', ['6']]
        Replace ['[*]'] with SM list.
        '''
        # for SMs
        sugar_indexes = [index for index, value in enumerate(template) if value == 'SUGARS']
        for i in sugar_indexes:
            template[i] = self.sugars
        full_indexes = [index for index, value in enumerate(template) if value == 'FULL_LIST']
        for j in full_indexes:
            template[j] = self.full_list
        # for sugars
        return template
