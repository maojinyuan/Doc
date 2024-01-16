import matplotlib.pyplot as plt
from matplotlib import font_manager
import os

# arial font_path = '/share/users/mjy/.fonts/arial.ttf'
# font_manager.fontManager.addfont(arial_font_path)
# plt.rcParams['font.family'] = ['Arial']
# plt.rcParams['font.sans-serif'] = ['Arial']

color_list = ['#8c564b', '#9467bd', '#1f77b4', '#2ca02c', '#e8572b', '#058b8c', '#800521', 'r']
fontsize=20
linesize=2
def make_figure(fontsize=20, linesize=2):
    plt.gcf().set_size_inches(8, 6)
    ax = plt.gca()

    ax.tick_params(labelsize=fontsize, width=linesize)
    ax.spines['bottom'].set_linewidth(linesize)
    ax.spines['left'].set_linewidth(linesize)
    ax.spines['top'].set_linewidth(linesize)
    ax.spines['right'].set_linewidth(linesize)
    ax.set_xlabel("$\overline{f}$", fontsize=fontsize)
    ax.set_ylabel("$\overline{x}_N/N$", fontsize=fontsize)
    # ax.set_title(os.path.basename(os.getcwd()), fontsize=fontsize)

    ax.tick_params(which='major', direction='in', length=8, width=linesize, labelsize=fontsize, right=True, top=True)
    ax.tick_params(which='minor', direction='in', length=4, width=linesize - 1, labelsize=fontsize, right=True, top=True)

    legend = ax.legend(fontsize=fontsize)
    legend.get_frame().set_linewidth(linesize)
    legend.get_frame().set_edgecolor('black')
    # plt.xlim(-0.1, 3.1)
    # plt.ylim(-0.1, 1.1)

def figure_show(save_dir=None, filename=None):
    save_dir = save_dir or os.getcwd()
    filename = filename or os.path.basename(save_dir)
    # save_path = os.path.join(str(save_dir), str(filename) + ".eps")
    save_path = os.path.join(str(save_dir), str(filename) + ".png")
    plt.savefig(save_path, dpi=600, bbox_inches='tight')


# def make_figure(fontsize=20, linesize=2):
    # plt.gcf().set_size_inches(8, 6)
    # ax = plt.gca()

    # ax.tick_params(labelsize=fontsize, width=linesize)
    # ax.spines['bottom'].set_linewidth(linesize)
    # ax.spines['left'].set_linewidth(linesize)
    # ax.spines['top'].set_linewidth(linesize)
    # ax.spines['right'].set_linewidth(linesize)

    # ax.set_xlabel("$\widetilde{f}$", fontsize=fontsize)
    # ax.set_ylabel("$\widetilde{x}_r$", fontsize=fontsize)
    # ax.set_title(os.path.basename(os.getcwd()), fontsize=fontsize)

    # ax.tick_params(which='major', direction='in', length=8, width=linesize, labelsize=fontsize, right=True, top=True)
    # ax.tick_params(which='minor', direction='in', length=4, width=linesize - 1, labelsize=fontsize, right=True, top=True)

    # legend = ax.legend(fontsize=fontsize)
    # legend.get_frame().set_linewidth(linesize)
    # legend.get_frame().set_edgecolor('black')

# # def make_figure(fontsize=20, linesize=2):
    # # plt.gcf().set_size_inches(8, 6)
    # # ax = plt.gca()

    # # ax.tick_params(labelsize=fontsize, width=linesize)
    # # ax.spines['bottom'].set_linewidth(linesize)
    # # ax.spines['left'].set_linewidth(linesize)
    # # ax.spines['top'].set_linewidth(linesize)
    # # ax.spines['right'].set_linewidth(linesize)
    # # ax.set_xlabel("$\widetilde{f}$", fontsize=fontsize)
    # # ax.set_ylabel("$\widetilde{x}_r$", fontsize=fontsize)
    # # ax.set_title(os.path.basename(os.getcwd()), fontsize=fontsize)

    # # ax.tick_params(which='major', direction='in', length=8, width=linesize, labelsize=fontsize, right=True, top=True)
    # # ax.tick_params(which='minor', direction='in', length=4, width=linesize - 1, labelsize=fontsize, right=True, top=True)

    # # legend = ax.legend(fontsize=fontsize)
    # # legend.get_frame().set_linewidth(linesize)
# #     legend.get_frame().set_edgecolor('black')

# def figure_show(save_dir=None, filename=None):
    # save_dir = save_dir or os.getcwd()
    # filename = filename or os.path.basename(save_dir)
    # save_path = os.path.join(str(save_dir), str(filename) + ".png")
    # plt.savefig(save_path, dpi=600, bbox_inches='tight')
