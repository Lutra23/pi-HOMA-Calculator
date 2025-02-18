import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import joblib
from typing import Dict, List, Tuple, Optional
import logging
from dataclasses import dataclass
import os

# 设置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# 理想键级值
IDEAL_BOND_ORDERS = {
    ('C', 'N'): 0.5129,
    ('N', 'C'): 0.5129,
    ('C', 'S'): 0.4497,
    ('S', 'C'): 0.4497,
    ('C', 'O'): 0.43615,
    ('O', 'C'): 0.43615,
    ('C', 'C'): 0.66,
}

@dataclass
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
        import shutil
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
        # 增加优化精度：
        # --opt tight: 使用更严格的优化标准
        # --acc 0.1: 提高SCF精度到0.1 meV
        # --cycles 1000: 增加最大优化步数
        # --maxiter 500: 增加最大SCF迭代次数
        # --etemp 300: 设置电子温度为300K以提高收敛性
        os.system('xtb xtbtopo.mol --opt tight --acc 0.1 --cycles 1000 --maxiter 500 --etemp 300 > xtbopt.log')
        
        # 读取优化后的结构
        if os.path.exists('xtbopt.mol'):
            optimized_mol = Chem.SDMolSupplier('xtbopt.mol')[0]
            # 清理临时文件
            os.remove('xtbtopo.mol')
            os.remove('xtbopt.mol')
            os.remove('xtbopt.log')
            if os.path.exists('charges'): os.remove('charges')
            if os.path.exists('wbo'): os.remove('wbo')
            if os.path.exists('xtbrestart'): os.remove('xtbrestart')
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

class MolecularFeatureExtractor:
    """分子特征提取器类"""
    
    # 电负性字典作为类变量
    ELECTRONEGATIVITY = {
        'H': 2.20, 'C': 2.55, 'N': 3.04, 
        'O': 3.44, 'F': 3.98, 'S': 2.58
    }
    
    def __init__(self):
        """初始化特征提取器"""
        self.logger = logging.getLogger(__name__)
        
    def generate_3d_mol(self, mol: Chem.Mol, mol_hash: str) -> Optional[Chem.Mol]:
        """生成分子的3D构象"""
        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            return mol
        except Exception as e:
            self.logger.error(f"生成3D构象失败 - Hash: {mol_hash}: {str(e)}")
            return None

    def get_atom_features(self, mol: Chem.Mol, atom_idx: int, prefix: str = '') -> Tuple[List[float], List[str]]:
        """提取原子特征"""
        atom = mol.GetAtomWithIdx(atom_idx)
        
        # 获取原子的环系统信息
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        rings_count = sum(1 for ring in atom_rings if atom_idx in ring)
        
        try:
            # 计算原子的Gasteiger电荷
            mol_copy = Chem.Mol(mol)
            AllChem.ComputeGasteigerCharges(mol_copy)
            gasteiger_charge = float(atom.GetProp('_GasteigerCharge'))
        except:
            gasteiger_charge = 0.0
        
        # 获取原子的邻居信息
        neighbors = atom.GetNeighbors()
        
        features = [
            # 基础原子属性
            float(atom.GetAtomicNum()),
            self.ELECTRONEGATIVITY.get(atom.GetSymbol(), 0.0),
            float(atom.GetDegree()),
            float(atom.GetHybridization()),
            float(atom.GetIsAromatic()),
            float(atom.GetFormalCharge()),
            float(atom.GetNumRadicalElectrons()),
            float(atom.IsInRing()),
            
            # 扩展原子属性
            float(atom.GetImplicitValence()),
            float(atom.GetExplicitValence()),
            float(atom.GetTotalValence()),
            float(atom.GetNumImplicitHs()),
            float(atom.GetNumExplicitHs()),
            float(atom.GetTotalNumHs()),
            float(atom.GetChiralTag()),
            
            # 电荷相关
            gasteiger_charge,
            atom.GetMass(),
            
            # 环系统相关
            float(rings_count),
        ]
        
        feature_names = [
            # 基础原子属性
            f'{prefix}_atomic_num',
            f'{prefix}_electronegativity',
            f'{prefix}_degree',
            f'{prefix}_hybridization',
            f'{prefix}_is_aromatic',
            f'{prefix}_formal_charge',
            f'{prefix}_num_radical_electrons',
            f'{prefix}_is_in_ring',
            
            # 扩展原子属性
            f'{prefix}_implicit_valence',
            f'{prefix}_explicit_valence',
            f'{prefix}_total_valence',
            f'{prefix}_num_implicit_hs',
            f'{prefix}_num_explicit_hs',
            f'{prefix}_total_num_hs',
            f'{prefix}_chiral_tag',
            
            # 电荷相关
            f'{prefix}_gasteiger_charge',
            f'{prefix}_mass',
            
            # 环系统相关
            f'{prefix}_rings_count',
        ]
        
        # 添加环的大小信息
        for ring_size in range(3, 8):
            in_ring = any(len(ring) == ring_size and atom_idx in ring for ring in atom_rings)
            features.append(float(in_ring))
            feature_names.append(f'{prefix}_in_ring{ring_size}')
        
        # 添加邻居信息
        num_neighbors = len(neighbors)
        num_c_neighbors = sum(1 for n in neighbors if n.GetSymbol() == 'C')
        num_n_neighbors = sum(1 for n in neighbors if n.GetSymbol() == 'N')
        num_o_neighbors = sum(1 for n in neighbors if n.GetSymbol() == 'O')
        num_f_neighbors = sum(1 for n in neighbors if n.GetSymbol() == 'F')
        num_s_neighbors = sum(1 for n in neighbors if n.GetSymbol() == 'S')
        num_aromatic_neighbors = sum(1 for n in neighbors if n.GetIsAromatic())
        
        features.extend([
            float(num_neighbors),
            float(num_c_neighbors),
            float(num_n_neighbors),
            float(num_o_neighbors),
            float(num_f_neighbors),
            float(num_s_neighbors),
            float(num_aromatic_neighbors),
            float(Chem.rdMolDescriptors.CalcHallKierAlpha(mol))  # Hall-Kier alpha值
        ])
        
        feature_names.extend([
            f'{prefix}_num_neighbors',
            f'{prefix}_num_c_neighbors',
            f'{prefix}_num_n_neighbors',
            f'{prefix}_num_o_neighbors',
            f'{prefix}_num_f_neighbors',
            f'{prefix}_num_s_neighbors',
            f'{prefix}_num_aromatic_neighbors',
            f'{prefix}_hall_kier_alpha'
        ])
        
        return features, feature_names

    def get_molecular_descriptors(self, mol: Chem.Mol) -> Tuple[List[float], List[str]]:
        """计算分子描述符"""
        try:
            from rdkit.Chem import Descriptors
            from rdkit.Chem.rdMolDescriptors import GetMACCSKeysFingerprint
            from rdkit.Chem.Descriptors3D import (
                NPR1, NPR2, PMI1, PMI2, PMI3, 
                Asphericity, Eccentricity, SpherocityIndex
            )
            
            # 基础描述符
            basic_descriptors = [
                (Descriptors.MolWt, 'mol_weight'),
                (Descriptors.MolLogP, 'logp'),
                (Descriptors.TPSA, 'tpsa'),
                (Descriptors.ExactMolWt, 'exact_mol_weight'),
                (Descriptors.FractionCSP3, 'fraction_csp3'),
                (Descriptors.HeavyAtomMolWt, 'heavy_atom_mol_weight'),
                (Descriptors.NumRotatableBonds, 'num_rotatable_bonds'),
                (Descriptors.NumHAcceptors, 'num_h_acceptors'),
                (Descriptors.NumHDonors, 'num_h_donors'),
                (Descriptors.NumHeteroatoms, 'num_heteroatoms'),
            ]
            
            # 形状描述符
            shape_descriptors = [
                (NPR1, 'npr1'),
                (NPR2, 'npr2'),
                (PMI1, 'pmi1'),
                (PMI2, 'pmi2'),
                (PMI3, 'pmi3'),
                (Asphericity, 'asphericity'),
                (Eccentricity, 'eccentricity'),
                (SpherocityIndex, 'spherocity_index')
            ]
            
            # 环系统描述符
            ring_descriptors = [
                (Descriptors.RingCount, 'ring_count'),
                (lambda m: sum(1 for x in range(m.GetNumAtoms()) if m.GetAtomWithIdx(x).GetIsAromatic()), 'aromatic_atoms_count'),
                (Descriptors.NumAromaticRings, 'num_aromatic_rings'),
            ]
            
            # 极性和电荷相关描述符
            AllChem.ComputeGasteigerCharges(mol)
            charge_descriptors = [
                (Descriptors.NumValenceElectrons, 'num_valence_electrons'),
                (lambda m: max(float(x.GetProp('_GasteigerCharge')) for x in m.GetAtoms() if '_GasteigerCharge' in x.GetPropsAsDict()), 'max_partial_charge'),
                (lambda m: min(float(x.GetProp('_GasteigerCharge')) for x in m.GetAtoms() if '_GasteigerCharge' in x.GetPropsAsDict()), 'min_partial_charge'),
            ]
            
            # 连接性描述符
            connectivity_descriptors = [
                (Descriptors.Chi0n, 'chi0n'),
                (Descriptors.Chi1n, 'chi1n'),
                (Descriptors.Chi2n, 'chi2n'),
                (Descriptors.Chi3n, 'chi3n'),
                (Descriptors.Chi4n, 'chi4n'),
                (Descriptors.Kappa1, 'kappa1'),
                (Descriptors.Kappa2, 'kappa2'),
                (Descriptors.Kappa3, 'kappa3'),
            ]
            
            # 计算所有描述符
            all_descriptors = (
                basic_descriptors + shape_descriptors + 
                ring_descriptors + charge_descriptors + 
                connectivity_descriptors
            )
            
            # 提取特征值和名称
            features = []
            names = []
            
            for descriptor_func, name in all_descriptors:
                try:
                    value = descriptor_func(mol)
                    features.append(float(value))
                    names.append(name)
                except Exception as e:
                    self.logger.warning(f"计算描述符 {name} 时出错: {str(e)}")
                    features.append(0.0)  # 使用默认值
                    names.append(name)
            
            return features, names
            
        except Exception as e:
            self.logger.error(f"计算分子描述符时出错: {str(e)}")
            return [], []

    def extract_features(self, mol: Chem.Mol, atom1_idx: int, atom2_idx: int) -> Tuple[List[float], List[str]]:
        """提取键和原子特征"""
        try:
            # 计算分子哈希值
            mol_hash = float(sum(ord(c) for c in Chem.MolToSmiles(mol)))
            
            # 获取键
            bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
            if bond is None:
                raise ValueError(f"原子之间不存在键: {atom1_idx}, {atom2_idx}")
            
            # 提取键特征
            bond_features = [
                float(bond.GetBondTypeAsDouble()),
                float(bond.GetIsConjugated()),
                float(bond.IsInRing()),
                float(bond.GetStereo()),
                abs(self.ELECTRONEGATIVITY.get(mol.GetAtomWithIdx(atom1_idx).GetSymbol(), 0.0) -
                    self.ELECTRONEGATIVITY.get(mol.GetAtomWithIdx(atom2_idx).GetSymbol(), 0.0))
            ]
            
            bond_feature_names = [
                'bond_type', 'is_conjugated', 'is_in_ring',
                'stereo', 'electronegativity_diff'
            ]
            
            # 提取原子特征
            atom1_features, atom1_names = self.get_atom_features(mol, atom1_idx, 'atom1')
            atom2_features, atom2_names = self.get_atom_features(mol, atom2_idx, 'atom2')
            
            # 计算原子间距离
            conf = mol.GetConformer()
            pos1 = conf.GetAtomPosition(atom1_idx)
            pos2 = conf.GetAtomPosition(atom2_idx)
            distance = float(np.linalg.norm([pos1.x - pos2.x, pos1.y - pos2.y, pos1.z - pos2.z]))
            
            # 提取分子特征和MACCS指纹
            mol_features, mol_names = self.get_molecular_descriptors(mol)
            maccs = [float(x) for x in list(Chem.rdMolDescriptors.GetMACCSKeysFingerprint(mol))]
            maccs_names = [f'MACCS_{i}' for i in range(len(maccs))]
            
            # 合并所有特征
            all_features = (
                [mol_hash] + bond_features + atom1_features + atom2_features +
                [distance] + mol_features + maccs
            )
            
            all_names = (
                ['mol_hash'] + bond_feature_names + atom1_names + atom2_names +
                ['distance'] + mol_names + maccs_names
            )
            
            return all_features, all_names
            
        except Exception as e:
            self.logger.error(f"提取特征时出错: {str(e)}")
            return None, None

class PiHOMACalculator:
    def __init__(self, model_path: str = 'best_gb_model.pkl', scaler_path: str = 'scaler.pkl'):
        """初始化π-HOMA计算器"""
        self.model = joblib.load(model_path)
        self.scaler = joblib.load(scaler_path)
        self.feature_extractor = MolecularFeatureExtractor()
        self.logger = logging.getLogger(__name__)

    def get_ring_bonds(self, mol: Chem.Mol) -> List[Tuple[int, int]]:
        """获取分子中所有环中的键"""
        ring_bonds = []
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            # 对环中的每个原子，检查它与环中下一个原子之间的键
            for i in range(len(ring)):
                atom1 = ring[i]
                atom2 = ring[(i + 1) % len(ring)]
                ring_bonds.append((atom1, atom2))
        return ring_bonds

    def predict_pi_bond_order(self, mol: Chem.Mol, atom1_idx: int, atom2_idx: int) -> Optional[float]:
        """预测两个原子之间的π键级"""
        try:
            # 生成3D构象
            mol_3d = self.feature_extractor.generate_3d_mol(mol, "")
            if mol_3d is None:
                return None

            # 提取特征
            features, _ = self.feature_extractor.extract_features(mol_3d, atom1_idx, atom2_idx)
            if features is None:
                return None

            # 标准化特征
            features_scaled = self.scaler.transform([features])

            # 预测π键级
            pi_bond_order = self.model.predict(features_scaled)[0]
            return pi_bond_order

        except Exception as e:
            self.logger.error(f"预测π键级时出错: {str(e)}")
            return None

    def calculate_pi_homa(self, mol: Chem.Mol, ring: List[int]) -> Optional[float]:
        """计算特定环的π-HOMA值"""
        try:
            n = len(ring)  # 环中键的数量
            normalized_deviations = []

            # 对环中的每个键计算π键级
            for i in range(n):
                atom1_idx = ring[i]
                atom2_idx = ring[(i + 1) % n]
                
                # 获取原子符号
                atom1_symbol = mol.GetAtomWithIdx(atom1_idx).GetSymbol()
                atom2_symbol = mol.GetAtomWithIdx(atom2_idx).GetSymbol()
                
                # 获取理想键级
                bond_type = (atom1_symbol, atom2_symbol)
                if bond_type not in IDEAL_BOND_ORDERS:
                    continue
                
                bo_ideal = IDEAL_BOND_ORDERS[bond_type]
                
                # 预测π键级
                pi_bond_order = self.predict_pi_bond_order(mol, atom1_idx, atom2_idx)
                if pi_bond_order is None:
                    continue
                
                # 计算归一化的偏差
                normalized_dev = (pi_bond_order - bo_ideal) ** 2 / (bo_ideal ** 2)
                normalized_deviations.append(normalized_dev)

            # 如果没有有效键，返回None
            if not normalized_deviations:
                return None

            # 计算π-HOMA
            # HOMA = 1 - (1/n) * Σ((BOi - BOideal)^2 / (BOideal^2))
            pi_homa = 1 - (sum(normalized_deviations) / len(normalized_deviations))
            return pi_homa

        except Exception as e:
            self.logger.error(f"计算π-HOMA时出错: {str(e)}")
            return None

    def calculate_all_rings_pi_homa(self, mol: Chem.Mol) -> Dict[int, float]:
        """计算分子中所有环的π-HOMA值
        Args:
            mol: RDKit分子对象
        Returns:
            Dict[int, float]: 环索引到π-HOMA值的映射
        """
        try:
            # 获取所有环
            ring_info = mol.GetRingInfo()
            rings = ring_info.AtomRings()

            # 计算每个环的π-HOMA
            results = {}
            for i, ring in enumerate(rings):
                pi_homa = self.calculate_pi_homa(mol, list(ring))
                if pi_homa is not None:
                    results[i] = pi_homa

            return results

        except Exception as e:
            self.logger.error(f"计算所有环的π-HOMA时出错: {str(e)}")
            return {}

def main():
    """主函数"""
    # 示例：计算分子的π-HOMA
    smiles = "n1ccccc1"  # 苯的SMILES表示
    
    try:
        # 创建分子对象（会自动进行xtb优化）
        molecule = MoleculeData(smiles)
            
        # 创建特征提取器
        extractor = MolecularFeatureExtractor()
        
        # 创建π-HOMA计算器
        calculator = PiHOMACalculator()
        
        # 计算所有环的π-HOMA
        results = calculator.calculate_all_rings_pi_homa(molecule.mol)
        
        # 打印结果
        print(f"分子SMILES: {smiles}")
        for ring_idx, pi_homa in results.items():
            print(f"环 {ring_idx + 1}:")
            # 获取环中的原子
            ring_atoms = list(molecule.mol.GetRingInfo().AtomRings()[ring_idx])
            print(f"  原子索引: {ring_atoms}")
            print(f"  π-HOMA值: {pi_homa:.4f}")
            print()
            
    except Exception as e:
        print(f"计算过程中出错: {str(e)}")

if __name__ == "__main__":
    main() 