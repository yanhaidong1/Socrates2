import pandas as pd
import os
import re
import argparse

# ========== 命令行参数解析 ==========
def parse_args():
    parser = argparse.ArgumentParser(description="处理共线区块文件")
    parser.add_argument('--file1', type=str, required=True, help="Pg__v__Pm.tsv 文件路径")
    parser.add_argument('--file2', type=str, required=True, help="Pm_vs_Pg.synHits.txt 文件路径")
    return parser.parse_args()

# ========== 步骤函数定义 ==========

def step1_filter_and_explode(file_path):
    print("步骤1：读取文件，过滤并拆分……")
    df = pd.read_csv(file_path, sep='\t', header=None)
    print(f"读取 {file_path}，共 {df.shape[0]} 行")
    df = df.drop(0, axis=1)
    filtered = df[
        (~df[1].astype(str).str.contains(',')) &
        (~df[2].astype(str).str.contains(','))
    ].copy()
    print(f"过滤后剩余 {filtered.shape[0]} 行")
    filtered[2] = filtered[2].astype(str).str.split(',')
    filtered[2] = filtered[2].apply(lambda x: [i.strip() for i in x])
    exploded = filtered.explode(2).reset_index(drop=True)

    out_path = file_path.replace('.tsv', '_filter_exploded.tsv')
    exploded.to_csv(out_path, sep='\t', index=False, header=False)
    print(f"步骤1完成，保存拆分文件：{out_path}\n")
    return out_path

def step2_filter_isAnchor(file_path):
    print("步骤2：筛选 isAnchor == True ……")
    base, ext = os.path.splitext(file_path)
    out_path = f"{base}_isAnchor_TRUE{ext}"
    df = pd.read_csv(file_path, sep='\t', header=0)
    anchor_df = df[df['isAnchor'] == True]
    anchor_df.to_csv(out_path, sep='\t', index=False)
    print(f"步骤2完成，保存结果：{out_path}\n")
    return out_path

def step3_grouping(exploded_path, is_anchor_path):
    print("步骤3：筛选连续区块并分组……")
    base3, ext3 = os.path.splitext(is_anchor_path)
    filtered_path = f"{base3}_filtered{ext3}"
    final_out = f"{base3}_filtered_final{ext3}"

    df_keys = pd.read_csv(exploded_path, sep='\t', usecols=[0], header=None)
    keys = set(df_keys.iloc[:, 0])
    df = pd.read_csv(is_anchor_path, sep='\t')
    filtered = df[df['id1'].isin(keys)].copy()

    filtered['ord1'] = pd.to_numeric(filtered['ord1'], errors='coerce')
    filtered['ord2'] = pd.to_numeric(filtered['ord2'], errors='coerce')
    filtered.dropna(subset=['ord1', 'ord2'], inplace=True)
    filtered.sort_values(by=['ord1', 'ord2'], inplace=True)
    filtered.reset_index(drop=True, inplace=True)

    filtered["duplicate"] = False
    for i in range(1, len(filtered) - 1):
        prev = filtered.loc[i - 1]
        curr = filtered.loc[i]
        next = filtered.loc[i + 1]
        prev_cont = abs(curr['ord1'] - prev['ord1']) == 1 and abs(curr['ord2'] - prev['ord2']) == 1
        next_cont = abs(next['ord1'] - curr['ord1']) == 1 and abs(next['ord2'] - curr['ord2']) == 1
        if prev_cont and next_cont:
            filtered.loc[i, "duplicate"] = True

    duplicated_rows = filtered[filtered["duplicate"]].copy()
    df_all = pd.concat([filtered, duplicated_rows], ignore_index=True)
    df_all.drop(columns="duplicate", inplace=True)
    df_all.sort_values(by=["ord1", "ord2"], inplace=True)
    df_all.reset_index(drop=True, inplace=True)

    new_groups = [0]
    group_id = 0
    for i in range(1, len(df_all)):
        prev = df_all.loc[i - 1]
        curr = df_all.loc[i]
        cont = abs(curr['ord1'] - prev['ord1']) == 1 and abs(curr['ord2'] - prev['ord2']) == 1
        if cont:
            new_groups.append(group_id)
        else:
            group_id += 1
            new_groups.append(group_id)
    df_all['group'] = new_groups

    continuous_groups = df_all.groupby('group').filter(lambda x: len(x) >= 2)
    continuous_groups.to_csv(final_out, sep='\t', index=False)
    print(f"步骤3完成，输出连续区块文件：{final_out}\n")
    return final_out

def step4_convert_fields(filtered_final_path):
    print("步骤4：转换字段，生成 list_name 等……")
    base, ext = os.path.splitext(filtered_final_path)
    out_path = base + '_converted.txt'

    df = pd.read_csv(filtered_final_path, sep='\t')
    df_out = pd.DataFrame()
    df_out['chr1'] = df['chr1']
    df_out['start1'] = df['start1']
    df_out['end1'] = df['end1']
    df_out['id1'] = df['id1']
    df_out['chr2'] = df['chr2']
    df_out['start2'] = df['start2']
    df_out['end2'] = df['end2']
    df_out['id2'] = df['id2']

    df_out['list_name'] = (
        df['genome1'].astype(str) + '_' +
        df['chr1'].astype(str) + '_' +
        df['start1'].astype(str) + '_' +
        df['end1'].astype(str) + '_' +
        df['genome2'].astype(str) + '_syntenic_region_' +
        df['group'].astype(str) + '_comparsion_' +
        df['blkID'].astype(str).str.replace(': ', '_', regex=False)
    )

    df_out['blkID'] = df['blkID']
    df_out['genome2'] = df['genome2']
    df_out['region_name'] = df['genome2'].astype(str) + '_syntenic_region_' + df['group'].astype(str)
    df_out['comparison_name'] = df['blkID'].astype(str).str.replace(': ', '_', regex=False)

    df_out.to_csv(out_path, sep='\t', index=False)
    print(f"步骤4完成，输出转换文件：{out_path}\n")
    return out_path

def step5_identify_gaps(converted_path):
    print("步骤5：识别间隙区段并保存 BED 文件……")
    df = pd.read_csv(converted_path, sep='\t')
    required_cols = {'chr2', 'start2', 'end2', 'region_name'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"输入文件缺少必要字段：{required_cols - set(df.columns)}")

    df_extracted = df[['chr2', 'start2', 'end2', 'region_name']]
    gap_rows = []

    for region, group in df_extracted.groupby('region_name'):
        group_sorted = group.sort_values('start2').reset_index(drop=True)
        for i in range(len(group_sorted) - 1):
            end_current = group_sorted.loc[i, 'end2']
            start_next = group_sorted.loc[i + 1, 'start2']
            if start_next > end_current:
                gap_rows.append({
                    'chr': group_sorted.loc[i, 'chr2'],
                    'start': end_current,
                    'end': start_next,
                    'region_name': region
                })

    gap_df = pd.DataFrame(gap_rows)
    if gap_df.empty:
        print("未检测到间隙区段，跳过生成 BED 文件。")
        return

    def get_prefix(region_name):
        match = re.match(r'([A-Za-z0-9]+_syntenic_region)', region_name)
        return match.group(1) if match else 'unknown_regions'

    output_dir = os.path.dirname(converted_path)
    for prefix, group in gap_df.groupby(gap_df['region_name'].apply(get_prefix)):
        bed_path = os.path.join(output_dir, f"{prefix}.bed")
        group[['chr', 'start', 'end', 'region_name']].to_csv(
            bed_path, sep='\t', index=False, header=False
        )
        print(f"保存 {len(group)} 条间隙区段到：{bed_path}")

    print("步骤5完成。\n")

# ========== 主程序入口 ==========

if __name__ == '__main__':
    args = parse_args()
    file_path_1 = args.file1
    file_path_2 = args.file2

    exploded_path = step1_filter_and_explode(file_path_1)
    is_anchor_path = step2_filter_isAnchor(file_path_2)
    filtered_final_path = step3_grouping(exploded_path, is_anchor_path)
    converted_path = step4_convert_fields(filtered_final_path)
    step5_identify_gaps(converted_path)
    print("所有步骤执行完毕！")
