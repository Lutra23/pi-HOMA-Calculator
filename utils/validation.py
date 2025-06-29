"""分子验证工具"""
import logging
import re
from typing import Tuple, Set
from rdkit import Chem
from rdkit.Chem import Descriptors
from config import LOG_CONFIG, SECURITY_CONFIG

# 配置日志
logging.basicConfig(**LOG_CONFIG)
logger = logging.getLogger(__name__)

def is_aromatic_ring(mol: Chem.Mol, ring_atoms: Tuple[int, ...]) -> bool:
    """
    判断环是否为芳香环
    
    Args:
        mol: RDKit分子对象
        ring_atoms: 环中的原子索引元组
    
    Returns:
        bool: 如果是芳香环返回True，否则返回False
    """
    try:
        # 检查环中的所有原子是否都是芳香的
        all_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms)
        if not all_aromatic:
            return False
            
        # 检查环中的所有键是否都是芳香的
        for i in range(len(ring_atoms)):
            atom1 = ring_atoms[i]
            atom2 = ring_atoms[(i + 1) % len(ring_atoms)]
            bond = mol.GetBondBetweenAtoms(atom1, atom2)
            if bond is None or not bond.GetIsAromatic():
                return False
        
        return True
    except Exception as e:
        logger.error(f"判断芳香环时出错: {str(e)}")
        return False

def validate_smiles_format(smiles: str) -> bool:
    """
    验证SMILES字符串格式
    
    Args:
        smiles: SMILES字符串
    
    Returns:
        bool: 格式是否有效
    """
    if not smiles or not isinstance(smiles, str):
        return False
    
    # 长度检查
    if len(smiles) > SECURITY_CONFIG['max_smiles_length']:
        return False
    
    # 基本字符检查
    allowed_chars = set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789()[]@+-=#$%./\\')
    if not all(c in allowed_chars for c in smiles):
        return False
    
    # 括号平衡检查
    if not _check_brackets_balance(smiles):
        return False
    
    return True

def _check_brackets_balance(smiles: str) -> bool:
    """检查括号是否平衡"""
    stack = []
    bracket_pairs = {'(': ')', '[': ']'}
    
    for char in smiles:
        if char in bracket_pairs:
            stack.append(char)
        elif char in bracket_pairs.values():
            if not stack:
                return False
            last_open = stack.pop()
            if bracket_pairs[last_open] != char:
                return False
    
    return len(stack) == 0

def validate_molecule_structure(mol: Chem.Mol) -> Tuple[bool, str]:
    """
    验证分子结构的合理性
    
    Args:
        mol: RDKit分子对象
    
    Returns:
        Tuple[bool, str]: (是否有效, 错误信息)
    """
    try:
        # 检查分子是否为空
        if mol is None:
            return False, "分子对象为空"
        
        # 检查原子数量
        num_atoms = mol.GetNumAtoms()
        if num_atoms == 0:
            return False, "分子不包含任何原子"
        
        if num_atoms > SECURITY_CONFIG['max_atoms']:
            return False, f"分子原子数量过多（{num_atoms} > {SECURITY_CONFIG['max_atoms']}）"
        
        # 检查原子类型
        allowed_atoms = SECURITY_CONFIG['allowed_atoms']
        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in allowed_atoms:
                return False, f"不支持的原子类型: {atom.GetSymbol()}"
        
        # 检查环数量
        ring_info = mol.GetRingInfo()
        num_rings = ring_info.NumRings()
        if num_rings > SECURITY_CONFIG['max_rings']:
            return False, f"环数量过多（{num_rings} > {SECURITY_CONFIG['max_rings']}）"
        
        # 检查分子量
        mol_weight = Descriptors.MolWt(mol)
        if mol_weight > 2000:  # 分子量上限
            return False, f"分子量过大: {mol_weight:.2f}"
        
        # 检查是否包含芳香环
        aromatic_rings = []
        for ring in ring_info.AtomRings():
            if is_aromatic_ring(mol, ring):
                aromatic_rings.append(ring)
        
        if not aromatic_rings:
            return False, "分子中未找到芳香环系"
        
        return True, ""
        
    except Exception as e:
        logger.error(f"验证分子结构时出错: {str(e)}")
        return False, f"分子结构验证失败: {str(e)}"

def validate_molecule_input(input_str: str) -> Tuple[Chem.Mol, str, str]:
    """
    验证分子输入并转换为RDKit分子对象
    
    Args:
        input_str: 输入的分子结构（SMILES或MOL格式）
    
    Returns:
        Tuple[Chem.Mol, str, str]: (RDKit分子对象, SMILES字符串, 输入格式)
        
    Raises:
        ValueError: 当输入无效时
    """
    if not input_str or not isinstance(input_str, str):
        raise ValueError("输入不能为空")
    
    input_str = input_str.strip()
    
    mol = None
    input_format = "unknown"
    smiles = ""
    
    try:
        # 检测输入格式
        if '\n' in input_str and ('M  END' in input_str or 'V2000' in input_str):
            # MOL格式
            mol = Chem.MolFromMolBlock(input_str)
            if mol is not None:
                input_format = "mol"
                smiles = Chem.MolToSmiles(mol)
        else:
            # SMILES格式
            if not validate_smiles_format(input_str):
                raise ValueError("SMILES格式无效")
            
            mol = Chem.MolFromSmiles(input_str)
            if mol is not None:
                input_format = "smiles"
                smiles = input_str
                
        if mol is None:
            raise ValueError("无法解析分子结构")
        
        # 验证分子结构
        is_valid, error_msg = validate_molecule_structure(mol)
        if not is_valid:
            raise ValueError(error_msg)
        
        # 标准化SMILES
        try:
            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            if canonical_smiles:
                smiles = canonical_smiles
        except Exception as e:
            logger.warning(f"SMILES标准化失败: {str(e)}")
        
        logger.info(f"分子验证成功 - 格式: {input_format}, SMILES: {smiles}")
        return mol, smiles, input_format
        
    except Exception as e:
        logger.error(f"分子验证失败: {str(e)}")
        raise ValueError(f"分子验证失败: {str(e)}")

def sanitize_input(input_str: str) -> str:
    """
    清理输入字符串
    
    Args:
        input_str: 原始输入字符串
    
    Returns:
        str: 清理后的字符串
    """
    if not input_str:
        return ""
    
    # 移除危险字符
    dangerous_chars = ['<', '>', '"', "'", '&', '\x00', '\r']
    cleaned = input_str
    for char in dangerous_chars:
        cleaned = cleaned.replace(char, '')
    
    # 限制长度
    max_length = SECURITY_CONFIG['max_smiles_length']
    if len(cleaned) > max_length:
        cleaned = cleaned[:max_length]
    
    return cleaned.strip()

def validate_calculation_parameters(params: dict) -> Tuple[bool, str]:
    """
    验证计算参数
    
    Args:
        params: 计算参数字典
    
    Returns:
        Tuple[bool, str]: (是否有效, 错误信息)
    """
    try:
        # 检查必需参数
        if 'smiles' not in params:
            return False, "缺少SMILES参数"
        
        smiles = params['smiles']
        if not smiles or not isinstance(smiles, str):
            return False, "SMILES参数无效"
        
        # 验证SMILES
        if not validate_smiles_format(smiles):
            return False, "SMILES格式无效"
        
        # 检查可选参数
        if 'options' in params:
            options = params['options']
            if not isinstance(options, dict):
                return False, "选项参数必须是字典类型"
            
            # 验证具体选项
            allowed_options = {
                'include_hydrogens', 'optimize_geometry', 
                'use_cache', 'timeout'
            }
            for key in options:
                if key not in allowed_options:
                    return False, f"不支持的选项: {key}"
        
        return True, ""
        
    except Exception as e:
        logger.error(f"验证计算参数时出错: {str(e)}")
        return False, f"参数验证失败: {str(e)}"

def get_molecule_info(mol: Chem.Mol) -> dict:
    """
    获取分子基本信息
    
    Args:
        mol: RDKit分子对象
    
    Returns:
        dict: 分子信息字典
    """
    try:
        info = {
            'num_atoms': mol.GetNumAtoms(),
            'num_bonds': mol.GetNumBonds(),
            'num_rings': mol.GetRingInfo().NumRings(),
            'molecular_weight': round(Descriptors.MolWt(mol), 2),
            'formula': Chem.rdMolDescriptors.CalcMolFormula(mol),
            'num_aromatic_rings': 0,
            'num_aromatic_atoms': 0
        }
        
        # 统计芳香环和芳香原子
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if is_aromatic_ring(mol, ring):
                info['num_aromatic_rings'] += 1
        
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                info['num_aromatic_atoms'] += 1
        
        return info
        
    except Exception as e:
        logger.error(f"获取分子信息时出错: {str(e)}")
        return {}