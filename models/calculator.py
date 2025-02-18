import logging
from typing import Dict, List, Optional, Tuple
import joblib
from rdkit import Chem
from .features import MolecularFeatureExtractor
from config import LOG_CONFIG, IDEAL_BOND_ORDERS, MODEL_PATH, SCALER_PATH

# 配置日志
logging.basicConfig(**LOG_CONFIG)
logger = logging.getLogger(__name__)

class PiHOMACalculator:
    """π-HOMA计算器类"""
    
    def __init__(self, model_path: str = MODEL_PATH, scaler_path: str = SCALER_PATH):
        """
        初始化π-HOMA计算器
        Args:
            model_path: 模型文件路径
            scaler_path: 标准化器文件路径
        """
        try:
            self.model = joblib.load(model_path)
            self.scaler = joblib.load(scaler_path)
        except Exception as e:
            logger.error(f"加载模型文件失败: {str(e)}")
            raise
            
        self.feature_extractor = MolecularFeatureExtractor()
        self.logger = logging.getLogger(__name__)

    def get_ring_bonds(self, mol: Chem.Mol) -> List[Tuple[int, int]]:
        """
        获取分子中所有环中的键
        Args:
            mol: RDKit分子对象
        Returns:
            List[Tuple[int, int]]: 环中键的原子索引对列表
        """
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
        """
        预测两个原子之间的π键级
        Args:
            mol: RDKit分子对象
            atom1_idx: 第一个原子的索引
            atom2_idx: 第二个原子的索引
        Returns:
            Optional[float]: π键级值，如果预测失败则返回None
        """
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
        """
        计算特定环的π-HOMA值
        Args:
            mol: RDKit分子对象
            ring: 环中原子的索引列表
        Returns:
            Optional[float]: π-HOMA值，如果计算失败则返回None
        """
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
        """
        计算分子中所有环的π-HOMA值
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