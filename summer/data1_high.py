import argparse
import os

# 定义判断read中"N"的数量是否过多的函数
def contains_too_many_N(read, n):
    return read.count("N") > n


# 定义判断read中低质量碱基比例是否过高的函数
def contains_too_many_low_quality_bases(read, quality, q, r):
    low_quality_bases = sum(1 for base_quality in quality if base_quality < q)
    return low_quality_bases / len(read) > r


def generate_high_quality_fastq(input_fname, output_fname, n, q, r):
    with open(input_fname) as fn:
        lines = fn.readlines()
    read_num = len(lines[::4])

    # 打开一个新的文件用于写入
    with open(output_fname, "w") as output_file:
        for i, line1, line2, line3, line4 in zip(range(read_num), lines[::4], lines[1::4], lines[2::4], lines[3::4]):
            quality = [ord(char) - 33 for char in line4.strip()]  # 转换质量分数
            # 如果read中"N"的数量过多或低质量碱基比例过高，则跳过该read
            if contains_too_many_N(line2, n) or contains_too_many_low_quality_bases(line2, quality, q, r):
                continue
            # 否则写入新的文件
            output_file.write(line1.strip() +'\n'+ line2.strip() +'\n'+ line3.strip() +'\n'+ line4.strip()+'\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a high quality FASTQ file by removing low quality reads")
    parser.add_argument("-i", "--input", help="Input FASTQ file", required=True)
    parser.add_argument("-o", "--output", help="Output FASTQ file", required=True)
    parser.add_argument("-n", type=int, help="Number of 'N' to remove a read", required=True)
    parser.add_argument("-q", type=int, help="Quality threshold to consider a base as low quality", required=True)
    parser.add_argument("-r", type=float, help="Ratio of low quality bases to remove a read", required=True)
    args = parser.parse_args()

    generate_high_quality_fastq(args.input, args.output, args.n, args.q, args.r)

# python data1_high.py -i data/data1_low.fq -o data/data1_high.fq -n 10 -q 20 -r 0.1