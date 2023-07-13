import matplotlib.pyplot as plt
import seaborn as sns
import random
import os
import data1_low
import data1_high


class fastqc:
    def __init__(self, filename):
        # 数据输入
         reads_dict = {}
         with open(filename) as fn:
            lines = fn.readlines()
            read_num = len(lines[::4])
            read_length = len(lines[1])
            for i, line1, line2, line4 in zip(range(read_num), lines[0::4], lines[1::4], lines[3::4]):
                if line1[0] == '@':
                    reads_dict[i + 1] = [line2.strip(), line4.strip()]
            self.read_dict = reads_dict
            self.read_num = read_num
            self.read_length = read_length
            self.filepath = filename
    # 可以加几个遍历函数，为后续增加功能提供便利

    def Get_GCcontent(self, str):
        return (str.count('C') + str.count('G')) / len(str)

    def ACTG_plot(self):
        data_dict = self.Get_ACTG()
        # 创建画布和子图
        fig, ax = plt.subplots()
        # 绘制曲线图
        for d in data_dict.values():
            ax.plot(range(1, self.read_length), d, linestyle='-')
        plt.legend(loc='upper right', frameon=False, labels=data_dict.keys())

        # 设置坐标轴标签和标题
        ax.set_xlabel('GC Content')
        ax.set_ylabel('Sequence')
        ax.set_title('Per Sequence GC Content')
        ax.set_xlim(0, self.read_length - 2)
        ax.set_ylim(0, 1)
        # 创建图例

        # 展示图形
        plt.show()

    def GC_plot(self):
        GC_dict = {}
        for key, read in self.read_dict.items():
            GC_dict[key] = self.Get_GCcontent(read[0])
        # 创建画布和子图
        fig, ax = plt.subplots()
        # 绘制曲线图
        sns.kdeplot(GC_dict.values(), linestyle='-', color='b')
        # 设置坐标轴标签和标题
        ax.set_ylabel('Sequence')
        ax.set_xlabel('GC Content')
        ax.set_title('Per Sequence GC Content')
        # 设置图例的位置和样式
        plt.legend = plt.legend(loc='upper right', frameon=False, labels=["GC Content"])
        # 展示图形
        plt.show()

    # 将ASCII编码转换为碱基质量
    def convert_to_quality(self, ascii_string):
        return [ord(char) - 33 for char in ascii_string]

    # 将ASCII编码转换为碱基质量
    def convert_to_ASCII(self, ascii_string):
        return [ord(char) + 33 for char in ascii_string]

    # 统计碱基质量分布
    def calculate_quality_distribution(self):
        quality_distribution = []
        for read in self.read_dict.values():
            qualities = self.convert_to_quality(read[1].strip())
            if len(quality_distribution) == 0:
                quality_distribution = [[q] for q in qualities]
            else:
                for j, q in enumerate(qualities):
                    quality_distribution[j].append(q)
        return quality_distribution

    # 绘制碱基质量分布图
    def plot_quality_distribution(self, quality_distribution):
        fig, ax = plt.subplots()
        pro_quality_distribution = []
        for read1, read2, read3 in zip(quality_distribution[::3], quality_distribution[1::3],
                                       quality_distribution[2::3]):
            pro_quality_distribution.append(read1 + read2 + read3)
        ax.boxplot(pro_quality_distribution, showfliers=False)
        ax.set_xlabel('Position')
        ax.set_ylabel('Base Quality')
        ax.set_title('Base Quality Distribution')
        # ax.legend()
        plt.show()

    # 替换某一行内容(某一行)
    def modify_file_line(self, file_path, line_number, new_content):
        # 读取文件内容
        with open(file_path, 'r') as file:
            lines = file.readlines()
        # 修改特定行的内容
        if line_number >= 1 and line_number <= len(lines):
            print(line_number)
            lines[line_number - 1] = str(new_content) + '\n'  # 添加换行符
        # 将修改后的内容写回文件
        with open(file_path, 'w') as file:
            file.writelines(lines)

    # 替换某一行内容(某一行); 可以加一个函数参数
    def modify_file_all(self, file_path):
        # 读取文件内容
        with open(file_path, 'r') as file:
            lines = file.readlines()
        # 统计所有 reads 在各位置上碱基质量分布
        for key, read in self.read_dict.items():
            lines[key*4-1] = str(self.convert_to_quality(read[1])) + '\n'
        # 将修改后的内容写回文件
        with open(file_path, 'w') as file:
            file.writelines(lines)

    def degrade_quality(self, p=0.05, k=15):
        for key, read in self.read_dict.items():
            # 以给定的概率 p 随机选择 reads
            if random.random() < p:
                # 从 read 中随机选择 k 个位置
                positions = random.sample(range(len(read[0])), k)

                # 将这 k 个位置上的碱基替换为字符"N"
                for pos in positions:
                    read[0] = read[0][:pos] + "N" + read[0][pos + 1:]
        return self.read_dict

    def Get_ACTG(self):
        # ACGT含量可视化
        ACTG_dict = {"A": [], "C": [], "T": [], "G": [], "N": []}
        readlines = [x[0] for x in self.read_dict.values()]
        for i in range(len(readlines[0])):
            loc = [x[i] for x in readlines]
            ACTG_dict["A"].append(loc.count("A") / self.read_num)
            ACTG_dict["C"].append(loc.count("C") / self.read_num)
            ACTG_dict["T"].append(loc.count("T") / self.read_num)
            ACTG_dict["G"].append(loc.count("G") / self.read_num)
            ACTG_dict["N"].append(loc.count("N") / self.read_num)
        return ACTG_dict

fname = f"data/data1_.fq"
fnameFor3 = f"data/data1.fq"
fname_low = f"data/data1_low.fq"
fname_high = f"data/data1_high.fq"

os.system("python data1_low.py -i data/data1_.fq -o data/data1_low.fq -p 0.05 -k 15")
os.system("python data1_high.py -i data/data1_low.fq -o data/data1_high.fq -n 10 -q 20 -r 0.1")

reads = fastqc(fname)
# 计算并可视化GC含量
# reads.GC_plot()
# ACGT含量可视化
# reads.ACTG_plot()
# 统计所有 reads 在各位置上碱基质量分布
# reads.plot_quality_distribution(reads.calculate_quality_distribution())
# reads.modify_file_all(fnameFor3)
# 统计低质量FASTQ文件
# reads_low = fastqc(fname_low)
# reads_low.ACTG_plot()
# 统计高质量FASTQ文件
reads_high = fastqc(fname_high)
reads_high.ACTG_plot()


