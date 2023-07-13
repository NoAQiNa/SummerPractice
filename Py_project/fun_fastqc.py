import matplotlib.pyplot as plt
import seaborn as sns
import random
import os
import data1_low
import data1_high

def DataInput(filename):
    # 数据输入
    reads_dict = {}
    with open(filename) as fn:
        lines = fn.readlines()
        read_num = len(lines[::4])
        read_length = len(lines[1])
        for i, line1, line2, line4 in zip(range(read_num), lines[0::4], lines[1::4], lines[3::4]):
            if line1[0] == '@':
                reads_dict[i + 1] = [line2.strip(), line4.strip()]
        return (reads_dict, read_num, read_length)

def Get_GCcontent(str):
    return (str.count('C') + str.count('G')) / len(str)

def plot(data_dict):
    # 创建画布和子图
    fig, ax = plt.subplots()

    # 绘制曲线图
    for d in data_dict.values():
        ax.plot(range(1,read_length), d, linestyle='-')
    plt.legend(loc='upper right', frameon=False, labels=data_dict.keys())

    # 设置坐标轴标签和标题
    ax.set_xlabel('GC Content')
    ax.set_ylabel('Sequence')
    ax.set_title('Per Sequence GC Content')
    ax.set_xlim(0, read_length-2)
    ax.set_ylim(0,1)
    # 创建图例

    # 展示图形
    plt.show()

def GC_plot(GC_dict):
    # 创建画布和子图
    fig, ax = plt.subplots()

    # 绘制曲线图
    sns.kdeplot(GC_dict.values(), linestyle='-', color='b')

    # 设置坐标轴标签和标题
    ax.set_ylabel('Sequence')
    ax.set_xlabel('GC Content')
    ax.set_title('Per Sequence GC Content')
    # 设置图例的位置和样式
    plt.legend = plt.legend(loc='upper right', frameon=False,labels=["GC Content"])
    # 展示图形
    plt.show()

# 将ASCII编码转换为碱基质量
def convert_to_quality(ascii_string):
    return [ord(char) - 33 for char in ascii_string]

# 统计碱基质量分布
def calculate_quality_distribution(fastq_lines):
    quality_distribution = []

    for i in range(3, len(fastq_lines), 4):
        qualities = convert_to_quality(fastq_lines[i].strip())
        if len(quality_distribution) == 0:
            quality_distribution = [[q] for q in qualities]
        else:
            for j, q in enumerate(qualities):
                quality_distribution[j].append(q)

    return quality_distribution

# 绘制碱基质量分布图
def plot_quality_distribution(quality_distribution):
    fig, ax = plt.subplots()
    pro_quality_distribution = []
    for read1,read2,read3 in zip(quality_distribution[::3],quality_distribution[1::3],quality_distribution[2::3]):
        pro_quality_distribution.append(read1+read2+read3)
    ax.boxplot(pro_quality_distribution, showfliers=False)
    ax.set_xlabel('Position')
    ax.set_ylabel('Base Quality')
    ax.set_title('Base Quality Distribution')
    # ax.legend()
    plt.show()

# 替换某一行内容
def modify_file_line(file_path, line_number, new_content):
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

def degrade_quality(reads_dict, p=0.05, k=15):
    for key, read in reads_dict.items():
        # 以给定的概率 p 随机选择 reads
        if random.random() < p:
            # 从 read 中随机选择 k 个位置
            positions = random.sample(range(len(read[0])), k)

            # 将这 k 个位置上的碱基替换为字符"N"
            for pos in positions:
                read[0] = read[0][:pos] + "N" + read[0][pos + 1:]
    return reads_dict

def Get_ACTG(reads_dict):
    # ACGT含量可视化
    ACTG_dict = {"A": [], "C": [], "T": [], "G": [], "N": []}
    readlines = [x[0] for x in reads_dict.values()]
    for i in range(len(readlines[0])):
        loc = [x[i] for x in readlines]
        ACTG_dict["A"].append(loc.count("A") / read_num)
        ACTG_dict["C"].append(loc.count("C") / read_num)
        ACTG_dict["T"].append(loc.count("T") / read_num)
        ACTG_dict["G"].append(loc.count("G") / read_num)
        ACTG_dict["N"].append(loc.count("N") / read_num)
    return ACTG_dict



fname = f"data/data1_.fq"
fnameFor3 = f"data/data1.fq"
fname_low = f"data/data1_low.fq"
fname_high = f"data/data1_high.fq"

os.system("python data1_low.py -i data/data1_.fq -o data/data1_low.fq -p 0.05 -k 15")
os.system("python data1_high.py -i data/data1_low.fq -o data/data1_high.fq -n 10 -q 20 -r 0.1")


(reads_dict, read_num, read_length) = DataInput(fname)
'''
# 计算并可视化GC含量
reads_GC_dict = {}
for key,read in reads_dict.items():
    reads_GC_dict[key] = Get_GCcontent(read[0])
GC_plot(reads_GC_dict)
'''
'''
# ACGT含量可视化
ACTG_dict = Get_ACTG(lines[1::4])
plot(ACTG_dict)
'''
'''
# 统计所有 reads 在各位置上碱基质量分布
for key, read in reads_dict.items():
    reads_dict[key][1] = convert_to_quality(read[1])
    modify_file_line(fnameFor3, key*4, reads_dict[key][1])

distribution = calculate_quality_distribution(lines)
plot_quality_distribution(distribution)
'''
'''
# 低质量文件
(reads_dict, read_num, read_length) = DataInput(fname_low)
ACTG_dict = Get_ACTG(reads_dict)
plot(ACTG_dict)
'''
'''
# 高质量文件ACGT含量可视化
(reads_dict, read_num, read_length) = DataInput(fname_high)
ACTG_dict = Get_ACTG(lines[1::4])
plot(ACTG_dict)
'''







