from rdkit import Chem
from rdkit.Chem import AllChem
import os
import logging
import shutil
from config import LOG_CONFIG, XTB_CONFIG

# 配置日志
logging.basicConfig(**LOG_CONFIG)
logger = logging.getLogger(__name__)

class MoleculeData:
    """分子数据的数据类"""
    def __init__(self, smiles: str):
        """
        初始化分子数据
        Args:
            smiles: SMILES字符串
        """
        self.logger = logging.getLogger(__name__)
        self.smiles = smiles
        
        # 创建分子并添加氢原子
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"无法从SMILES创建分子: {smiles}")
        
        # 添加氢原子
        self.mol = Chem.AddHs(mol)
        
        # 生成3D构象并优化
        AllChem.EmbedMolecule(self.mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(self.mol)
        
        # 检查是否安装了xtb
        if shutil.which('xtb') is not None:
            # 保存初始xyz文件
            self._save_xyz("initial.xyz")
            # 使用xtb优化
            self.mol = self.optimize_with_xtb(self.mol)
        else:
            self.logger.warning("未找到xtb程序，将使用RDKit的MMFF94优化结果")

    def _save_xyz(self, filename: str):
        """保存分子为xyz格式"""
        with open(filename, 'w') as f:
            # 写入原子数
            f.write(f"{self.mol.GetNumAtoms()}\n")
            # 写入注释行
            f.write(f"Generated from SMILES: {self.smiles}\n")
            # 写入坐标
            conf = self.mol.GetConformer()
            for i in range(self.mol.GetNumAtoms()):
                atom = self.mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

    def optimize_with_xtb(self, mol):
        """使用xtb优化分子构型"""
        # 将分子写入临时文件
        writer = Chem.SDWriter("xtbtopo.mol")
        writer.write(mol)
        writer.close()
        
        # 运行xtb优化
        xtb_command = (
            f'xtb xtbtopo.mol '
            f'--opt {XTB_CONFIG["opt_level"]} '
            f'--acc {XTB_CONFIG["accuracy"]} '
            f'--cycles {XTB_CONFIG["max_cycles"]} '
            f'--maxiter {XTB_CONFIG["max_iterations"]} '
            f'--etemp {XTB_CONFIG["electronic_temp"]} '
            '> xtbopt.log'
        )
        os.system(xtb_command)
        
        # 读取优化后的结构
        if os.path.exists('xtbopt.mol'):
            optimized_mol = Chem.SDMolSupplier('xtbopt.mol')[0]
            # 清理临时文件
            temp_files = [
                'xtbtopo.mol', 'xtbopt.mol', 'xtbopt.log',
                'charges', 'wbo', 'xtbrestart'
            ]
            for file in temp_files:
                if os.path.exists(file):
                    os.remove(file)
            return optimized_mol
        else:
            raise Exception("xtb优化失败")

    def _load_xyz_coordinates(self, xyz_file: str):
        """从xyz文件加载原子3D坐标"""
        try:
            # 读取xyz文件
            with open(xyz_file, 'r') as f:
                lines = f.readlines()
                
            # 第一行是原子数
            n_atoms = int(lines[0].strip())
            
            # 跳过注释行
            coords = []
            atoms = []
            for line in lines[2:n_atoms+2]:
                atom, x, y, z = line.strip().split()
                atoms.append(atom)
                coords.append([float(x), float(y), float(z)])
                
            # 创建新的构象对象
            conf = Chem.Conformer(self.mol.GetNumAtoms())
            
            # 设置坐标
            for i, (x, y, z) in enumerate(coords):
                conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(x, y, z))
                
            # 将构象添加到分子中
            self.mol.RemoveAllConformers()  # 移除所有现有构象
            self.mol.AddConformer(conf)     # 添加新构象
                
        except Exception as e:
            raise ValueError(f"加载xyz坐标时出错: {str(e)}")