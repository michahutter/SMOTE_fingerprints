import scipy.stats as stats
import matplotlib.pyplot as plt

from source.code import utils


def create_evaluation_output(group1, group2, name1, name2, out_file, sig_level=0.05):
    """ Used to do a visual and statistical evaluation of the results found in two compare files.

    For the visual evaluation, the average similarity scores found in two compare files are transformed in two separate
    histograms. Having the possible scores on the x_axis and the amount with which those appear on the y_axis. The
    values are rounded to one decimal, such that the histograms are more compact.
    For the statistical evaluation a Mann-Whitney-U test is performed on the specific values found in the two compared
    files. The p-value of the test will be stored in a file as well as the test_statistic U.

    Parameters
    ----------
    group1 : list of list
        Contains the data found in the first compare file.
    group2 : list of list
        Contain the data found in the second compare file.
    name1 : str
        Name of the first file to be used.
    name2 :
        Name of the second file to be used.
    out_file : str
        Name to be used for the out file.
    sig_level : float
        Significance level that should be used to determine at which p-value the null-hypothesis is to be rejected.

    Returns
    -------
    int
        Return 1 if the null hypothesis was rejected, 1 if it could not be rejected, None if errors occurred.
    """

    out_file = utils.check_extension_output(out_file, '.asc')
    if out_file is None:
        return None

    group1_val = [float(i[1]) for i in group1]
    group2_val = [float(i[1]) for i in group2]

    res = stats.mannwhitneyu(group1_val, group2_val, alternative='greater')

    if res.pvalue < sig_level:
        body = 'The calculated p-value: ' + str(res.pvalue) + ' is smaller than the chosen cut off: ' + \
               str(sig_level) + '\n\nH0-Hypothesis is rejected\nTest statistic: ' + str(res.statistic)
        res = 1

    else:
        body = 'The calculated p-value: ' + str(res.pvalue) + ' is bigger than the chosen cut off: ' + \
               str(sig_level) + '\n\nH0-Hypothesis cannot be rejected\nTest statistic: ' + str(res.statistic)
        res = 0

    write_evaluation_file(out_file, name1, name2, body)

    # Creation of bar charts for compare files
    outfile1 = out_file[:-4] + '_hist1.png'
    outfile2 = out_file[:-4] + '_hist2.png'
    barchart(group1, outfile1)
    barchart(group2, outfile2)

    return res


def write_evaluation_file(out_file, name1, name2, body):
    """ Helper function, takes care of writing the output file for create_evaluation_output()

    Parameters
    ----------
    out_file : str
        output file specified in create_evaluation_output().
    name1 : str
        Name of first compare file specified in create_evaluation_output().
    name2 : str
        Name of second compare file specified in create_evaluation_output().
    body : str
        Text to be used in the output file specified in create_evaluation_output().
    """

    with open(out_file, 'w+') as f:
        # Header of out file
        f.write('Evaluation of similarity scores between ' + name1 + ' and ' + name2 + ':\n\n\n')
        f.write('H0: The similarity scores are equal between the two groups\n')
        f.write('H1: The similarity scores are larger for: ' + name1 + '\n\n')
        f.write(body)
        f.close()


def barchart(data, out_file):
    """ Helper function, used to create the bar charts for create_evaluation_output()

    Parameters
    ----------
    data : list of list
        Every list contains a list with an identifier and a similarity score.
    out_file : str
        Name for generated output specified in create_evaluation_output().
    """

    if out_file is None:
        return None

    bins = list()
    x = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

    scores = list()
    for i in range(len(data)):
        score = round(float(data[i][1]), 1)
        scores.append(score)

    bins.append(round((x[0] - 0.05), 2))

    for i in range(len(x)):
        bins.append(round((x[i] + 0.05), 2))

    x_label = 'Similarity'
    y_label = 'Amount'

    counts, edges, bars = plt.hist(scores, color='royalblue', bins=bins, rwidth=0.9)
    plt.title('Distribution of similarities')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xticks(x)

    try:
        plt.bar_label(bars, labels=[f'{d:.1E}' if int(d) > 1000000 else f'{d:.0f}' for d in bars.datavalues])
    except:
        print('WARNING: Matplotlib version to low, using simpler barchart formatting.')

    plt.savefig(out_file)
    plt.close()


def formatting(x):
    """
    Helper function for creation of bar charts.

    Parameters
    ----------
    x : list of str

    Returns
    -------
    list of str
        Formatted string.
    """
    for i in range(len(x)):
        if int(x[i]) > 1000000:
            x = format(x, '%.2E')

    return x
