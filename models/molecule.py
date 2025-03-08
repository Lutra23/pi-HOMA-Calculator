from rdkit import Chem
from rdkit.Chem import AllChem
import os
import logging
import shutil
import hashlib
import time
import tempfile
from pathlib import Path
from typing import Optional, Dict, Any, Tuple
from functools import lru_cache
from config import LOG_CONFIG, XTB_CONFIG, TEMP_DIR, CACHE_DIR, cached

# 配置日志
logging.basicConfig(**LOG_CONFIG)
logger = logging.getLogger(__name__)

# 分子缓存字典和最大缓存大小
_molecule_cache: Dict[str, Any] = {}
MAX_CACHE_SIZE = 50  # 最大缓存分子数量

class MoleculeData:
    """分子数据的数据类，用于处理分子的3D结构生成和优化"""
    
    def __init__(self, smiles: str):
        """
        初始化分子数据
        
        Args:
            smiles: SMILES字符串表示的分子结构
            
        Raises:
            ValueError: 当SMILES字符串无法转换为分子结构时
            RuntimeError: 当分子结构生成或优化失败时
        """
        self.logger = logging.getLogger(__name__)
        self.smiles = smiles
        self.mol_hash = self._generate_hash(smiles)
        
        # 尝试从缓存获取分子
        if self.mol_hash in _molecule_cache:
            self.logger.info(f"从缓存获取分子: {smiles}")
            self.mol = _molecule_cache[self.mol_hash]
            return
            
        # 如果缓存过大，清理最早的条目
        if len(_molecule_cache) >= MAX_CACHE_SIZE:
            # 简单策略：删除第一个键
            try:
                first_key = next(iter(_molecule_cache))
                del _molecule_cache[first_key]
                self.logger.debug(f"缓存已满，删除条目: {first_key}")
            except Exception as e:
                self.logger.warning(f"清理缓存时出错: {str(e)}")
            
        # 尝试从文件缓存获取分子
        cache_file = os.path.join(CACHE_DIR, f"{self.mol_hash}.mol")
        if os.path.exists(cache_file):
            try:
                self.logger.info(f"从文件缓存加载分子: {smiles}")
                self.mol = Chem.SDMolSupplier(cache_file)[0]
                if self.mol is not None:
                    _molecule_cache[self.mol_hash] = self.mol
                    return
            except Exception as e:
                self.logger.warning(f"从缓存加载分子失败: {str(e)}，将重新生成")
        
        # 创建分子并添加氢原子
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"无法从SMILES创建分子: {smiles}")
            
            # 添加氢原子
            self.mol = Chem.AddHs(mol)
            
            # 生成3D构象并优化
            self._generate_3d_structure()
            
            # 缓存分子
            _molecule_cache[self.mol_hash] = self.mol
            
            # 保存到文件缓存
            self._save_to_cache(cache_file)
            
        except Exception as e:
            self.logger.error(f"初始化分子数据时出错: {str(e)}")
            raise RuntimeError(f"初始化分子数据时出错: {str(e)}")
    
    def _generate_hash(self, smiles: str) -> str:
        """为SMILES字符串生成唯一哈希值"""
        return hashlib.md5(smiles.encode()).hexdigest()
    
    def _save_to_cache(self, cache_file: str) -> None:
        """将分子保存到文件缓存"""
        try:
            writer = Chem.SDWriter(cache_file)
            writer.write(self.mol)
            writer.close()
            self.logger.debug(f"分子已保存到缓存: {cache_file}")
        except Exception as e:
            self.logger.warning(f"保存分子到缓存失败: {str(e)}")
    
    def _generate_3d_structure(self) -> None:
        """生成分子的3D结构并优化"""
        try:
            # 生成初始3D构象
            AllChem.EmbedMolecule(self.mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(self.mol)
            
            # 如果启用了XTB优化且XTB可用，则使用XTB进一步优化
            if XTB_CONFIG['enabled'] and shutil.which('xtb') is not None:
                self.logger.info("使用XTB优化分子构型")
                try:
                    # 创建临时目录进行优化，避免文件冲突
                    with tempfile.TemporaryDirectory(dir=TEMP_DIR) as temp_dir:
                        # 保存初始xyz文件
                        xyz_file = os.path.join(temp_dir, "initial.xyz")
                        self._save_xyz(xyz_file)
                        # 使用xtb优化
                        self.mol = self._optimize_with_xtb(self.mol, temp_dir)
                except Exception as e:
                    self.logger.warning(f"XTB优化失败，回退到MMFF优化: {str(e)}")
                    # 如果XTB优化失败，确保分子有有效的3D构象
                    if not self.mol.GetNumConformers():
                        self.logger.info("重新生成MMFF优化的3D构象")
                        AllChem.EmbedMolecule(self.mol, randomSeed=42)
                        AllChem.MMFFOptimizeMolecule(self.mol)
            else:
                if not XTB_CONFIG['enabled']:
                    self.logger.info("XTB优化已禁用，使用RDKit的MMFF94优化结果")
                else:
                    self.logger.warning("未找到xtb程序，使用RDKit的MMFF94优化结果")
        except Exception as e:
            self.logger.error(f"生成3D结构时出错: {str(e)}")
            raise RuntimeError(f"生成3D结构时出错: {str(e)}")

    def _save_xyz(self, filename: str) -> None:
        """保存分子为xyz格式
        
        Args:
            filename: 要保存的文件路径
        """
        try:
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
            self.logger.debug(f"分子已保存为XYZ格式: {filename}")
        except Exception as e:
            self.logger.error(f"保存XYZ文件时出错: {str(e)}")
            raise IOError(f"保存XYZ文件时出错: {str(e)}")

    def _optimize_with_xtb(self, mol: Chem.Mol, temp_dir: str) -> Chem.Mol:
        """使用xtb优化分子构型
        
        Args:
            mol: 要优化的分子
            temp_dir: 临时目录路径
            
        Returns:
            优化后的分子
            
        Raises:
            RuntimeError: 当XTB优化失败时
        """
        # 文件路径
        input_file = os.path.join(temp_dir, "xtbtopo.mol")
        output_file = os.path.join(temp_dir, "xtbopt.mol")
        log_file = os.path.join(temp_dir, "xtbopt.log")
        
        try:
            # 将分子写入临时文件
            writer = Chem.SDWriter(input_file)
            writer.write(mol)
            writer.close()
            
            # 构建XTB命令
            xtb_command = (
                f'cd {temp_dir} && xtb xtbtopo.mol '
                f'--opt {XTB_CONFIG["opt_level"]} '
                f'--acc {XTB_CONFIG["accuracy"]} '
                f'--cycles {XTB_CONFIG["max_cycles"]} '
                f'--maxiter {XTB_CONFIG["max_iterations"]} '
                f'--etemp {XTB_CONFIG["electronic_temp"]} '
                f'> {log_file}'
            )
            
            # 运行XTB优化
            start_time = time.time()
            exit_code = os.system(xtb_command)
            end_time = time.time()
            
            self.logger.info(f"XTB优化完成，耗时: {end_time - start_time:.2f}秒")
            
            # 检查优化是否成功
            if exit_code != 0:
                with open(log_file, 'r') as f:
                    log_content = f.read()
                self.logger.error(f"XTB优化失败，退出码: {exit_code}\n日志内容: {log_content}")
                raise RuntimeError(f"XTB优化失败，退出码: {exit_code}")
            
            # 读取优化后的结构
            if os.path.exists(output_file):
                optimized_mol = Chem.SDMolSupplier(output_file)[0]
                if optimized_mol is None:
                    raise RuntimeError("无法读取XTB优化后的分子结构")
                return optimized_mol
            else:
                raise RuntimeError("XTB优化后的文件不存在")
                
        except Exception as e:
            self.logger.error(f"XTB优化过程中出错: {str(e)}")
            raise RuntimeError(f"XTB优化过程中出错: {str(e)}")

    def _load_xyz_coordinates(self, xyz_file: str) -> None:
        """从xyz文件加载原子3D坐标
        
        Args:
            xyz_file: XYZ文件路径
            
        Raises:
            ValueError: 当XYZ文件格式错误或无法读取时
        """
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
            
            self.logger.debug(f"已从XYZ文件加载坐标: {xyz_file}")
                
        except Exception as e:
            self.logger.error(f"加载XYZ坐标时出错: {str(e)}")
            raise ValueError(f"加载XYZ坐标时出错: {str(e)}")

def clear_molecule_cache():
    """清理分子缓存"""
    global _molecule_cache
    cache_size = len(_molecule_cache)
    _molecule_cache.clear()
    logger.info(f"已清理{cache_size}个分子缓存条目")

@cached
def get_molecule(smiles: str) -> Chem.Mol:
    """获取分子对象，带有缓存功能
    
    Args:
        smiles: SMILES字符串
        
    Returns:
        RDKit分子对象
        
    Raises:
        RuntimeError: 当分子创建失败时
    """
    try:
        molecule_data = MoleculeData(smiles)
        return molecule_data.mol
    except Exception as e:
        logger.error(f"获取分子对象失败: {str(e)}")
        raise RuntimeError(f"获取分子对象失败: {str(e)}")

def is_valid_smiles(smiles: str) -> bool:
    """检查SMILES字符串是否有效
    
    Args:
        smiles: 要检查的SMILES字符串
        
    Returns:
        如果SMILES有效则返回true，否则返回false
    """
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None