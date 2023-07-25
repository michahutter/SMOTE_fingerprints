import argparse
from os.path import exists

from source.code import fingerprints
from source.code import combine
from source.code import distances
from source.code.ranks import *
from source.code import statistics
from source.code import matching
from source.code import compare
from source.code import evaluation
from source.code import update_database
from source.code import update_pharma_fps
from source.code import update_join
from source.code import update_fps
from source.code import create_db


def main():
    """ Function for handling the input from the command line.

    Takes the input flags and calls the corresponding functions to fulfill the user requests.
    """

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-h', '--help', nargs='*')
    parser.add_argument('-f', '--fingerprints', nargs='*', type=str)
    parser.add_argument('-d', '--distances', nargs='*', type=str)
    parser.add_argument('-r', '--ranks', nargs='*', type=str)
    parser.add_argument('-c', '--combine', nargs='*')
    parser.add_argument('-s', '--statistic', nargs='*')
    parser.add_argument('-m', '--matching', nargs='*', type=str)
    parser.add_argument('-p', '--compare', nargs='*', type=str)
    parser.add_argument('-e', '--evaluation', nargs='*', type=str)
    parser.add_argument('-i', '--initialize', nargs='*')

    args = parser.parse_args()

    # If no flag was set, print general help text.
    if (args.fingerprints is None) and (args.distances is None) and (args.ranks is None) and \
            (args.combine is None) and (args.help is None) and (args.statistic is None) and \
            (args.matching is None) and (args.compare is None) and (args.evaluation is None) and \
            (args.initialize is None):

        print_help()

    # Fingerprint flag parsing.
    elif args.fingerprints is not None:
        a = args.fingerprints
        success = False

        if len(a) == 0:
            success = True
            print_help_file('source/help_texts/fingerprints_help.txt')
            return

        elif len(a) == 3:
            success = True
            res = fingerprints.create_fingerprints_output(a[0], a[1], a[2])
            if res is None:
                return

        elif len(a) == 4:
            if (a[3].lower() == 'Torsion'.lower()) or (a[3].lower() == 'Pairs'.lower()) or \
                    (a[3].lower() == 'Pharmacophore'.lower()):
                success = True
                res = fingerprints.create_fingerprints_output(a[0], a[1], a[2] + ' ' + a[3])

                if res is None:
                    return

        elif len(a) == 5:
            if a[3].lower() == 'Pharmacophore'.lower():
                success = True
                res = fingerprints.create_fingerprints_output(a[0], a[1], a[2] + ' ' + a[3], a[4])
                if res is None:
                    return
        else:
            print('\nERROR: Flags -f, --fingerprints require 3/4 or 5 arguments.\n')
            return

        if success:
            print('\nGeneration of fingerprints successful.\n')
        else:
            print('\nERROR: Invalid input arguments.\n')


    # Distance flag parsing.
    elif args.distances is not None:
        a = args.distances

        if len(a) == 0:
            print_help_file('source/help_texts/distances_text.txt')
            return

        elif len(a) == 3:
            actives = fingerprints.read_byte_fingerprint_file(a[0])
            inactives = fingerprints.read_byte_fingerprint_file(a[1])
            if actives is None or inactives is None:
                return

            res = distances.create_distances_output(actives[3], inactives[3], actives[1], inactives[1], a[2])
            if res is not None:
                print('\nCalculation of distances successful.\n')

        else:
            print('\nERROR: Flags -d, --distances require 3 arguments\n')
            return

    # Ranking flag parsing.
    elif args.ranks is not None:
        a = args.ranks

        if len(a) == 0:
            print_help_file('source/help_texts/ranks_help.txt')
            return

        data = distances.read_distances(a[0])
        if data is None:
            return

        if len(a) == 2:
            ranking = ranking_manhatten(data)
        elif (len(a) == 3) & (a[2].lower() == 'invert'.lower()):
            ranking = ranking_manhatten(data, False)
        else:
            print('\nERROR: Flags -r, --rank require 2 or 3 arguments\n')
            return

        create_ranking_output(ranking, a[1])
        print('\nRanking of compounds successful.')

    # Combine flag parsing.
    elif args.combine is not None:
        a = args.combine

        if len(a) == 0:
            print_help_file('source/help_texts/combine_help.txt')
            return

        elif len(a) == 5 or len(a) == 6:

            ranks = read_ranking_output(a[0])
            fps = fingerprints.read_byte_fingerprint_file(a[1])
            fp_type = fingerprints.read_byte_header_type(a[1])

            if ranks is None or fps is None:
                return

            if len(a) == 6:
                # Is given argument a valid float number?
                if a[5].replace('.', '', 1).isdigit():
                    combine.create_combine_output(ranks, fps, fp_type, a[2], a[3], int(a[4]), float(a[5]))
                else:
                    print('\nERROR: The value: ' + a[6] + ' is not a valid significance level')
                    return

            if len(a) == 5:
                combine.create_combine_output(ranks, fps, fp_type, a[2], a[3], int(a[4]))

        else:
            print('\nERROR: Flags -c, --combine require 4/5 or 6 arguments\n')
            return

        print('\nCombination of fingerprints successful.')

    # Statistic flag parsing.
    elif args.statistic is not None:
        a = args.statistic

        if len(a) == 0:
            print_help_file('source/help_texts/statistic_help.txt')
            return

        elif len(a) == 2:
            data_raw = combine.read_byte_fingerprint_file(a[0])
            if data_raw is None:
                return

            statistics.create_statistics_output(data_raw, a[1])

        else:
            print('\nERROR: Flags -s, --statistic require 2 arguments\n')
            return

        print('\nCreation of statistic successful.')

    # Matching flag parsing.
    elif args.matching is not None:
        a = args.matching
        if len(a) == 0:
            print_help_file('source/help_texts/matching_help.txt')
            return

        [data, names] = combine.read_combine(a[0])
        com_fp = fingerprints.regenerate_fingerprint(data[1], data[0])

        # Compare with SureChEMBL database
        if len(a) == 2:
            res = matching.get_matches(com_fp, data[0], a[1], names)
        elif len(a) == 3:
            if a[2].lower() == 'true':
                res = matching.get_matches(com_fp, data[0], a[1], names, True)
            elif a[2].lower() == 'false':
                res = matching.get_matches(com_fp, data[0], a[1], names, False)
            elif exists(a[2]):
                res = matching.get_matches_external_base(com_fp, data[0], a[1], a[2])
            else:
                print('\nERROR: Invalid arguments given.')
                return

        else:
            print('\nERROR: Flags -m, --matching require 2 or 3 arguments.')
            return

        if res:
            print('\nMatching of combined fingerprint to database successful.')

    # Compare flag parsing.
    elif args.compare is not None:
        a = args.compare
        if len(a) == 0:
            print_help_file('source/help_texts/compare_help.txt')
            return
        if len(a) == 3:
            res = matching.read_matching_file(a[0])
            if res is None:
                return

            [fp_type, schembls] = res
            compare_set = fingerprints.read_byte_fingerprint_file(a[1])
            if compare_set is None:
                return

            compare.compare(fp_type, schembls, compare_set[3], compare_set[1], a[2])
            print('\nComparing of compounds to used SureChEMBL entries successful.')

        # In case of external database
        elif len(a) == 4:
            res = matching.read_matching_file(a[0])
            if res is None:
                return

            [_, identifiers] = res
            compare_set = fingerprints.read_byte_fingerprint_file(a[1])
            if compare_set is None:
                return

            d_base = fingerprints.read_byte_fingerprint_file(a[3])
            compare.compare_external(identifiers, compare_set[3], compare_set[1], a[2], d_base[3], d_base[1])
            print('\nComparing of compounds to used database successful.')

        else:
            print('\nERROR: Flags -p, --compare require 3 arguments.')
            return

    # Evaluation flag parsing.
    elif args.evaluation is not None:
        a = args.evaluation
        if len(a) == 0:
            print_help_file('source/help_texts/evaluation_help.txt')
            return

        elif len(a) == 2:
            data = matching.read_all_scores(a[0])
            if data is None:
                return

            evaluation.barchart(data, a[1])
            print("\nCreation of all_scores bar chart successful.")

        elif len(a) == 3:
            group1 = compare.read_compare_file(a[0])
            group2 = compare.read_compare_file(a[1])
            if group1 is None or group2 is None:
                return

            res = evaluation.create_evaluation_output(group1, group2, a[0], a[1], a[2])

            if res is not None:
                print('\nEvaluation of similarities determined with compare command successful.')

        elif len(a) == 4:
            group1 = compare.read_compare_file(a[0])
            group2 = compare.read_compare_file(a[1])
            sig_level = float(a[3])

            # Is given argument a valid float number?
            if a[3].replace('.', '', 1).isdigit():

                res = evaluation.create_evaluation_output(group1, group2, a[0], a[1], a[2], sig_level)
                if res is not None:
                    print('\nEvaluation of similarities determined with compare command successful.')

            else:
                print('\nERROR: The value: ' + a[3] + ' is not a valid significance level')
                return

        else:
            print('\nERROR: Flags -e, --evaluation require 2 or 3/4 arguments.')
            return

    # Update flag parsing.
    elif args.initialize is not None:

        print('This process is used to recreate the contents of the local MYSQL database used for the thesis. It\n'
              'takes several days to complete this task for the drug-like compounds. It is not advised to generate\n'
              'all fingerprints for all compounds, since that will take several weeks to complete. Moreover, the\n'
              'matching processes only considers the filtered compounds by default, which was also done for the\n'
              'results in the thesis. Since only compounds with drug-like features are really of interest.\n\n')

        print('To continue with only generating the fingerprints of filtered compounds press "y", to continue with\n'
              'generating all press "n", to exit this command without any effects press "e".\n\n')

        ans = input()

        if ans.lower() == 'y':
            create_db.create_db()
            update_database.update()
            update_join.update()
            update_fps.update()
            update_pharma_fps.update(True)

        elif ans.lower() == 'n':
            create_db.create_db()
            update_database.update()
            update_join.update()
            update_fps.update()
            update_pharma_fps.update(False)

        elif ans.lower() == 'e':
            print('\nUpdating command exited.')
            return

        print('\nUpdating of database completed.')

    # Help flag parsing.
    elif args.help is not None:
        print_help()


# Help texts needed to be displayed.
def print_help():
    """ Function used for printing the general help text to the console.

    This function is used when the help flag was set or the program was started without any input arguments.
    The text is structured in a way that help texts for specific commands can be easily changed if needed. This
    can be done by specifying the contents of those help texts below.
    """

    print("USAGE: main.py [-h] [--help]\n"
          "               [-f] [--fingerprints [in_file out_file fp_type]]\n"
          "               [-d] [--distances    [fp_file1 fp_file2 out_file]]\n"
          "               [-r] [--ranks        [d_file out_file (invert)]]\n"
          "               [-c] [--combine      [r_file fp_file out_file mode amount (threshold)]]\n"
          "               [-s] [--statistic    [fp_file out_file]]\n"
          "               [-m] [--matching     [comb_file out_file (True/False)]]\n"
          "               [-m] [--matching     [comb_file out_file dataset_file]]\n"
          "               [-p] [--compare      [match_file fp_file out_file]]\n"
          "               [-p] [--compare      [match_file fp_file out_file dataset_fp_file]]\n"
          "               [-e] [--evaluation   [comp_file1 comp_file2 out_file (sig_level)]]\n"
          "               [-e] [--evaluation   [all_scores_file1 out_file]]\n"
          "               [-i] [--initialize]\n\n"

          + description_text + "\n"
          "OPTIONS:\n   -h, --help\n                        "
          + help_text +
          "\n   -f, --fingerprints [in_file out_file fingerprint_type (optional)]\n                        "
          + fp_text_general +
          "\n   -d, --distances [in_file1 in_file2 out_file]\n                        "
          + distances_text_general +
          "\n   -r, --ranks [in_file out_file (optional:invert)]\n                        "
          + ranks_text_general +
          "\n   -c, --combine [rank_file fp_file out_file mode amount (optionals:threshold)]\n"
          + "                        "
          + combine_text_general +
          "\n   -s, --statistic [fp_file out_file]\n                        "
          + statistic_text_general +
          "\n   -m, --matching [combine_file out_file (optional:true/false)]                        " +
          "\n   -m, --matching [combine_file out_file dataset_file]\n                        "
          + matching_text_general1 +
          "\n   -p, --compare [combine_file out_file dataset_file]]\n                        "
          + compare_text_general +
          "\n   -e, --evaluation [compare_file1 compare_file2 out_file (optional:float 0 - 1)]"
          + "                        " +
          "\n   -e, --evaluation [all_scores_file1 out_file]\n"
          + "                        "
          + evaluation_text_general1 +
          "\n   -i, --initialize\n                        "
          + initialize_text_general)


def print_help_file(file):
    """ Function used for printing command specific help texts to the console.

    Parameters
    ----------
    file : str
        Name of the file containing help text.
    """

    print('\n')
    with open(file, 'r') as file:
        text = file.read()
    print(text)


description_text = "This program is used to generate a combined fingerprint from a group of actives\n" \
                   "and inactives towards a specific target and compare it to a database, to find \n" \
                   "new possible substances with an effect on that target. The chosen compounds re-\n" \
                   "present the ones with properties found the least in the corresponding group of \n" \
                   "inactives for that target. To gain more information about specific commands run \n" \
                   "them without any input arguments.\n"

fp_text_general = "This command is used to create an output file contain-\n" \
                  "                        ing the original SMILES representation of the compound,\n" \
                  "                        the name of the compound, the size of the serialized ge-\n" \
                  "                        nerated fingerprints and the serialized byte representa-\n" \
                  "                        tion of the created fingerprints.\n"

distances_text_general = "Takes two files generated by the fingerprints command\n" \
                         "                        and calculates the manhatten distances of the two\n" \
                         "                        groups. An output file containing the names and dis-\n" \
                         "                        tances will be created. Actives must be used as the\n" \
                         "                        first argument.\n"\

ranks_text_general = "Takes a file generated by the distances command and\n" \
                     "                        summarizes the manhatten distances in a sum. Results\n" \
                     "                        will be grouped by the names in the first column. Out-\n" \
                     "                        put file will be in descending order by default, to \n" \
                     "                        get ascending order add 'invert' as extra last input\n" \
                     "                        argument.\n" \

combine_text_general = "This command is used to create a combined fingerprint\n" \
                       "                        file. It requires a ranking file generated with the\n" \
                       "                        ranks command and a fingerprint file generated with\n" \
                       "                        the fingerprints command. It is possible to control\n" \
                       "                        what kind of fingerprints and how many should be con-\n" \
                       "                        sidered, when creating the combined fingerprint. Op-\n" \
                       "                        tionally threshold can also be changed. Also a statis-\n" \
                       "                        tics file with fitting barchart is created.\n" \

statistic_text_general = "This command can be used to generate a statistics\n" \
                         "                        file, in the same style of the ones generated in the\n"\
                         "                        combine command. It contains information about set\n" \
                         "                        bits in a list of given fingerprints. It requires a\n"\
                         "                        fingerprint file generated with the the fingerprints\n" \
                         "                        command.\n" \

matching_text_general1 = "This command is used to match a combined fingerprint\n" \
                        "                        to the SureChEMBL data. File containing the combined\n" \
                        "                        fingerprint has to be generated by the combine com-\n" \
                        "                        mand. A file containing the top 10 most similar com-\n" \
                        "                        pounds will be created. Showing PubChem SID, Sure-\n" \
                        "                        ChEMBL ID, similarity score and patent IDs. Filter-\n"\
                        "                        ing of compounds can be turned off.\n"\
                        "                        Alternatively, it can be used to match a combined\n" \
                        "                        fingerprint to an additionally provided database. This\n" \
                        "                        database has to be in the form of a file, that follows\n" \
                        "                        the format of a fingerprints file or a compounds file.\n" \
                        "                        However, no additional info can be gathered from the\n" \
                        "                        PubChem services and no filtering possible.\n" \



compare_text_general = "This command is used to compare the initial compounds\n" \
                       "                        to the top 10 matches found via the matching command.\n" \
                       "                        The average similarity using the Tanimoto coefficient,\n" \
                       "                        between the given fingerprints and the found top 10\n" \
                       "                        matches will be calculated. A file containing the\n" \
                       "                        average of the resulting similarities to all given\n" \
                       "                        fingerprints in descending order will be generated.\n" \
                       "                        Average of all scores will be given at the end of the\n" \
                       "                        file.\n" \


evaluation_text_general1 = "This command is used to perform an evaluation of\n" \
                          "                        different kinds Depending on the input. Given two\n"\
                          "                        compare files, a Mann-Whitney-U test will be performed\n" \
                          "                        on them. The significance level can be given as an ad-\n" \
                          "                        ditional argument, otherwise 0.05 will be used. More- \n" \
                          "                        over for each compare file, a barchart showing the dis-\n" \
                          "                        tribution of similarity scores will be created. If one\n" \
                          "                        all_scores file is provided, a bar char representing the\n" \
                          "                        distribution will be created.\n"


initialize_text_general = "This command creates the local SureChEMBL database that\n" \
                      "                        is used for matching. This will take some time so only\n" \
                      "                        use it when no usage of the other commands is needed for\n" \
                      "                        the rest of the day.\n"

help_text = "Show this help message and exit\n"

if __name__ == "__main__":
    main()
    # print(combine.create_fingerprints_output.__doc__)
