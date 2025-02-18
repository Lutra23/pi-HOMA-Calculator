from typing import List, Tuple, Optional
import numpy as np
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import GetMACCSKeysFingerprint
from rdkit.Chem.Descriptors3D import (
    NPR1, NPR2, PMI1, PMI2, PMI3,
    Asphericity, Eccentricity, SpherocityIndex
)
from rdkit.Chem import Descriptors
from config import LOG_CONFIG, ELECTRONEGATIVITY

# 配置日志
logging.basicConfig(**LOG_CONFIG)
logger = logging.getLogger(__name__)

class MolecularFeatureExtractor:
    """分子特征提取器类"""
    
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
            ELECTRONEGATIVITY.get(atom.GetSymbol(), 0.0),
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
                abs(ELECTRONEGATIVITY.get(mol.GetAtomWithIdx(atom1_idx).GetSymbol(), 0.0) -
                    ELECTRONEGATIVITY.get(mol.GetAtomWithIdx(atom2_idx).GetSymbol(), 0.0))
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
            maccs = [float(x) for x in list(GetMACCSKeysFingerprint(mol))]
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