import argparse
import random
import os



def degrade_quality(input_fname, output_fname, p, k):
    with open(input_fname) as fn:
        lines = fn.readlines()
    read_num = len(lines[::4])
    read_length = len(lines[1])

    # 打开一个新的文件用于写入
    with open(output_fname, "w") as output_file:
        for i, line1, line2, line3, line4 in zip(range(read_num), lines[::4], lines[1::4], lines[2::4], lines[3::4]):
            # 以给定的概率 p 随机选择 reads
            if random.random() < p:
                # 从 read 中随机选择 k 个位置
                positions = random.sample(range(len(line2)), k)
                # 将这 k 个位置上的碱基替换为字符"N"
                for pos in positions:
                    line2 = line2[:pos] + "N" + line2[pos + 1:]
            # 写入修改后的read
            output_file.write(line1.strip() +'\n'+ line2.strip() +'\n'+ line3.strip() +'\n'+ line4.strip()+'\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a low quality FASTQ file by randomly replacing bases with 'N'")
    parser.add_argument("-i", "--input", help="Input FASTQ file", required=True)
    parser.add_argument("-o", "--output", help="Output FASTQ file", required=True)
    parser.add_argument("-p", type=float, help="Probability to select a read", required=True)
    parser.add_argument("-k", type=int, help="Number of positions to replace with 'N'", required=True)
    args = parser.parse_args()

    degrade_quality(args.input, args.output, args.p, args.k)

# python data1_low.py -i data/data1_.fq -o data/data1_low.fq -p 0.05 -k 15