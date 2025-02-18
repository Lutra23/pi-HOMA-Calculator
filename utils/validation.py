"""分子验证工具"""
import logging
from rdkit import Chem
from config import LOG_CONFIG

# 配置日志
logging.basicConfig(**LOG_CONFIG)
logger = logging.getLogger(__name__)

def is_aromatic_ring(mol, ring_atoms):
    """
    判断环是否为芳香环
    
    Args:
        mol: RDKit分子对象
        ring_atoms: 环中的原子索引列表
    
    Returns:
        bool: 如果是芳香环返回True，否则返回False
    """
    # 检查环中的所有原子是否都是芳香的
    all_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms)
    if not all_aromatic:
        return False
        
    # 检查环中的所有键是否都是芳香的
    for i in range(len(ring_atoms)):
        atom1 = ring_atoms[i]
        atom2 = ring_atoms[(i + 1) % len(ring_atoms)]
        bond = mol.GetBondBetweenAtoms(atom1, atom2)
        if not bond.GetIsAromatic():
            return False
    return True

def validate_molecule_input(input_str: str) -> tuple[Chem.Mol, str, str]:
    """
    验证分子输入并转换为RDKit分子对象
    
    Args:
        input_str: 输入的分子结构（SMILES或MOL格式）
    
    Returns:
        tuple: (RDKit分子对象, SMILES字符串, 输入格式)
        
    Raises:
        ValueError: 当输入无效时
    """
    mol = None
    input_format = "unknown"
    smiles = ""
    
    try:
        if '\n' in input_str:  # 可能是MOL格式
            mol = Chem.MolFromMolBlock(input_str)
            if mol is not None:
                input_format = "mol"
                smiles = Chem.MolToSmiles(mol)
        
        if mol is None:  # 如果不是MOL格式或转换失败，尝试SMILES格式
            mol = Chem.MolFromSmiles(input_str)
            if mol is not None:
                input_format = "smiles"
                smiles = input_str
                
        if mol is None:
            raise ValueError("无效的分子结构")
            
        # 验证分子中是否存在芳香环
        ring_info = mol.GetRingInfo()
        rings = list(ring_info.AtomRings())
        has_aromatic_ring = any(is_aromatic_ring(mol, ring) for ring in rings)
        
        if not has_aromatic_ring:
            raise ValueError("分子中未找到芳香环系")
            
        return mol, smiles, input_format
        
    except Exception as e:
        logger.error(f"分子验证失败: {str(e)}")
        raise ValueError(f"分子验证失败: {str(e)}")