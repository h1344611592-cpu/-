# 导入需要的库和工具函数
import argparse
import pandas as pd
from tools import (
    appendPrediction,
    append_structure_images,
    chinese_to_english_punctuation,
    clean_chemical_name,
    default_process,
    insert_images_to_excel,
    replace_greek_letters,
    
)
from tqdm import tqdm


# 启用进度条
tqdm.pandas()


# 兼容命令行参数
def parse_arguments():
    """
    解析命令行输入参数，供命令行调用使用。
    """
    parser = argparse.ArgumentParser(description="化合物列表处理")
    parser.add_argument("-i", "--input", type=str, default="data/testdata.xlsx", help="化合物列表输入路径")
    parser.add_argument("-o", "--output", type=str, default="data/testresult.xlsx", help="结果文件输出路径")
    parser.add_argument("-d", "--image_dir", type=str, default="test_imgs", help="图片保存目录")
    parser.add_argument("-n", "--names_col", type=str, default="英文名", help="英文名列名称")
    parser.add_argument("-s", "--step", type=str, choices=["step1", "step2", "all"], default="all", help="选择需要执行的步骤"),
    parser.add_argument("-hh", "--has_header", type=str, choices=["yes", "no"], default="yes",
                        help="输入文件是否包含列名（yes 或 no，默认 yes）")
    parser.add_argument("-ii", "--insert_images", action="store_true", help="是否向excel中插入结构图")
    parser.add_argument("-ss", "--suppress_state", action="store_true", help="是否输出中间提示信息")
    return parser.parse_args()


# 定义一般流程
def main_run(input_path="data/testdata.xlsx",
             output_path="data/result.xlsx",
             image_dir="test_imgs",
             names_col="英文名",
             suppress_state=False,
             insert_images=False,
             has_header = "yes",
             step = "all"):
    """
    主流程：输入英文名文件，输出带有CID、SMILES、预测结果和结构式图片的Excel文件。
    
    参数:
        - input_path: str，输入文件的路径（Excel）
        - output_path: str，输出文件的路径（Excel）
        - image_dir: str，保存结构图片的文件夹
        - names_col: str，输入文件中存放英文名的列名
        - insert_images: bool，是否在Excel中插入结构图
        - has_header: str，是否有列名（"yes" 或 "no"）
    """
    # 文件读取和英文名格式化
    header = 0 if args.has_header == "yes" else None
    df_input = pd.read_excel(input_path, header=header)
    
    if header is None:
        col_count = df_input.shape[1]
        if col_count == 1:
            df_input.columns = [names_col]
        else:
            df_input.columns = [names_col] + [f"Unnamed_{i}" for i in range(1, col_count)]
    
    # Step 1：CID和SMILES抓取
    if step in ("step1", "all"):
        print("正在执行 Step 1: 获取CID和SMILES...")
        
        df_input[names_col] = df_input[names_col].apply(chinese_to_english_punctuation)
        df_input[names_col] = df_input[names_col].astype(str).apply(clean_chemical_name)
        df_input[names_col] = df_input[names_col].apply(replace_greek_letters)
        print("英文名已格式化")
        
        # step 1：cid和smiles抓取
        df_input[["CID", "SMILES"]] = df_input[names_col].progress_apply(
            lambda name: default_process(name, keep_file=False, suppress_state=suppress_state)
        ).apply(pd.Series)
        
        # 重复抓取未找到的cid和smiles，排除网络问题
        retry_mask = df_input["CID"].isin(["NotFound", "MultipleCIDs", None]) | df_input["SMILES"].isin(["NoStructure", None])
        
        if retry_mask.sum() > 0:
            print(f"准备重试 {retry_mask.sum()} 条未成功的记录。")
            df_retry = df_input.loc[retry_mask].copy()
            df_retry[["CID", "SMILES"]] = df_retry[names_col].progress_apply(
                lambda name: default_process(name, keep_file=False, suppress_state=suppress_state)
            ).apply(pd.Series)
            
            df_input.update(df_retry)  # 更新数据框
            
        # 如果只执行 Step 1，提前保存结果
        if step == "step1":
            df_input = append_structure_images(df_input, smiles_col="SMILES", image_dir=image_dir)
            df_input.to_excel(output_path, index=False)
            print("Step 1 完成！输出文件：", output_path)
            return
        
    # Step 2：跨环境模型预测和结构图生成
    if step in ("step2", "all"):
        print("正在执行 Step 2: 模型预测和结构图生成...")
        # 如果从 Step 2 开始，需确保输入数据包含 CID 和 SMILES
        if "CID" not in df_input.columns or "SMILES" not in df_input.columns:
            raise ValueError("执行 Step 2 需要输入数据包含 'CID' 和 'SMILES' 列！")
        
        # step 2：跨环境模型预测
        df_result = appendPrediction(df_input, suppress_state=suppress_state)
        # 生成结构图
        df_result = append_structure_images(df_result, smiles_col="SMILES", image_dir=image_dir)
        # 写出结果文件
        df_result.to_excel(output_path, index=False)
        # 向excel中插图
        if insert_images:
            insert_images_to_excel(output_path, image_column="StructureImage")
        
        print("处理完成！输出文件：", output_path)  # 提示完成
        
# 实际运行
if __name__ == "__main__":
    args = parse_arguments()
    
    main_run(
        input_path=args.input,
        output_path=args.output,
        image_dir=args.image_dir,
        names_col=args.names_col,
        suppress_state=args.suppress_state,
        insert_images = args.insert_images,
        has_header=args.has_header,
        step=args.step
    )



