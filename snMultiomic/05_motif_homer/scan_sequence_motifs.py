#!/usr/bin/env python
"""
扫描一段序列中的motif

使用方法:
    python scan_sequence_motifs.py --sequence "ATCGATCGATCG" --motif_db path/to/motif.jaspar
    或者
    python scan_sequence_motifs.py --sequence_file sequence.fasta --motif_db path/to/motif.jaspar
"""

import argparse
import tempfile
from pathlib import Path

import pandas as pd
from scprinter import motifs
from scprinter import genome as genome_module


def scan_sequence_for_motifs(
    sequence: str,
    motif_db_path: str,
    genome_fa_path: str = None,
    pvalue: float = 5e-5,
    clean: bool = True,
    strand: bool = True,
    verbose: bool = True,
):
    """
    扫描序列中的motif
    
    Parameters
    ----------
    sequence : str
        要扫描的DNA序列（只包含ACGT字符）
    motif_db_path : str
        Motif数据库文件路径（JASPAR格式）
    genome_fa_path : str, optional
        基因组fasta文件路径。如果为None，会创建一个临时文件
    pvalue : float, optional
        P-value阈值，默认5e-5
    clean : bool, optional
        是否合并重叠的motif hit，默认True
    strand : bool, optional
        是否考虑链信息，默认True
    verbose : bool, optional
        是否显示进度，默认True
    
    Returns
    -------
    pd.DataFrame
        包含motif扫描结果的DataFrame，列包括：
        - chrom: 染色体名
        - start: 序列起始位置（相对于输入序列）
        - end: 序列结束位置（相对于输入序列）
        - tf_name: 转录因子名称
        - score: motif匹配分数
        - strand: 链信息（+或-）
        - motif_start: motif在序列中的起始位置
        - motif_end: motif在序列中的结束位置
    """
    # 验证序列
    sequence = sequence.upper().strip()
    valid_bases = set("ACGTN")
    if not all(c in valid_bases for c in sequence):
        invalid_chars = set(sequence) - valid_bases
        raise ValueError(f"序列包含无效字符: {invalid_chars}")
    
    # 如果序列包含N，提示用户
    if "N" in sequence:
        print("警告: 序列中包含N字符，可能会影响motif扫描结果")
    
    # 创建临时fasta文件（如果未提供基因组文件）
    temp_fa = None
    if genome_fa_path is None:
        temp_fa = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
        temp_fa.write(f">temp_sequence\n{sequence}\n")
        temp_fa.close()
        genome_fa_path = temp_fa.name
    
    try:
        # 初始化Motifs类
        motif_scanner = motifs.Motifs(
            ref_path_motif=motif_db_path,
            ref_path_fa=genome_fa_path,
            bg="even",  # 均匀背景频率
            pvalue=pvalue,
            n_jobs=1,  # 单个序列不需要多进程
        )
        
        # 准备scanner
        motif_scanner.prep_scanner(pvalue=pvalue)
        
        # 扫描motif（使用临时染色体名称）
        peaks = [["temp_sequence", 0, len(sequence)]]
        results = motif_scanner.scan_motif(
            peaks,
            clean=clean,
            concat=True,
            verbose=verbose,
            strand=strand,
        )
        
        # 转换为DataFrame
        if len(results) == 0:
            print("未发现任何motif匹配")
            return pd.DataFrame(columns=[
                "chrom", "start", "end", "index", "tf_name", 
                "score", "strand", "motif_start", "motif_end"
            ])
        
        df = pd.DataFrame(
            results,
            columns=[
                "chrom", "start", "end", "index", "tf_name",
                "score", "strand", "motif_start", "motif_end"
            ]
        )
        
        # 将链信息转换为更易读的格式
        df["strand"] = df["strand"].map({1: "+", -1: "-", "*": "*"})
        
        return df
        
    finally:
        # 清理临时文件
        if temp_fa is not None and Path(temp_fa.name).exists():
            Path(temp_fa.name).unlink()


def main():
    parser = argparse.ArgumentParser(
        description="扫描序列中的motif",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
    # 直接扫描序列
    python scan_sequence_motifs.py --sequence "ATCGATCGATCG" --motif_db motifs.jaspar
    
    # 从文件读取序列
    python scan_sequence_motifs.py --sequence_file sequence.fasta --motif_db motifs.jaspar
    
    # 指定p-value阈值
    python scan_sequence_motifs.py --sequence "ATCGATCGATCG" --motif_db motifs.jaspar --pvalue 1e-4
    
    # 保存结果到CSV
    python scan_sequence_motifs.py --sequence "ATCGATCGATCG" --motif_db motifs.jaspar -o results.csv
        """
    )
    
    # 输入选项（互斥）
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--sequence", "-s",
        type=str,
        help="要扫描的DNA序列"
    )
    input_group.add_argument(
        "--sequence_file", "-f",
        type=str,
        help="包含序列的FASTA文件路径"
    )
    
    # 必需参数
    parser.add_argument(
        "--motif_db", "-m",
        type=str,
        required=True,
        help="Motif数据库文件路径（JASPAR格式）"
    )
    
    # 可选参数
    parser.add_argument(
        "--genome_fa", "-g",
        type=str,
        default=None,
        help="基因组fasta文件路径（如果提供，可以用于计算背景频率）"
    )
    
    parser.add_argument(
        "--pvalue", "-p",
        type=float,
        default=5e-5,
        help="P-value阈值（默认: 5e-5）"
    )
    
    parser.add_argument(
        "--no-clean",
        action="store_true",
        help="不合并重叠的motif hit"
    )
    
    parser.add_argument(
        "--no-strand",
        action="store_true",
        help="不考虑链信息"
    )
    
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None,
        help="输出CSV文件路径（如果不指定，将打印到标准输出）"
    )
    
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="不显示进度信息"
    )
    
    args = parser.parse_args()
    
    # 读取序列
    if args.sequence:
        sequence = args.sequence
    else:
        with open(args.sequence_file, 'r') as f:
            lines = f.readlines()
            # 简单解析FASTA格式
            sequence = ""
            for line in lines:
                line = line.strip()
                if line and not line.startswith(">"):
                    sequence += line
            if not sequence:
                parser.error(f"无法从文件 {args.sequence_file} 中读取序列")
    
    # 扫描motif
    try:
        df = scan_sequence_for_motifs(
            sequence=sequence,
            motif_db_path=args.motif_db,
            genome_fa_path=args.genome_fa,
            pvalue=args.pvalue,
            clean=not args.no_clean,
            strand=not args.no_strand,
            verbose=not args.quiet,
        )
        
        # 输出结果
        if args.output:
            df.to_csv(args.output, index=False)
            print(f"\n结果已保存到: {args.output}")
        else:
            print("\n=== Motif扫描结果 ===")
            print(df.to_string(index=False))
            print(f"\n总共发现 {len(df)} 个motif匹配")
            if len(df) > 0:
                print(f"\n发现的转录因子: {', '.join(sorted(df['tf_name'].unique()))}")
        
    except Exception as e:
        parser.error(f"错误: {e}")


if __name__ == "__main__":
    main()

