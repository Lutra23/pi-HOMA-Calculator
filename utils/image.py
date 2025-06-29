"""分子图像生成工具"""
import logging
from typing import List, Tuple, Optional
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor, rdMolDraw2D
from rdkit import Chem
from config import LOG_CONFIG, MOLECULE_IMAGE_CONFIG

# 配置日志
logging.basicConfig(**LOG_CONFIG)
logger = logging.getLogger(__name__)

def generate_molecule_image(mol: Chem.Mol, rings: List[Tuple[int, ...]], 
                          width: Optional[int] = None, 
                          height: Optional[int] = None) -> bytes:
    """
    生成带有环高亮的分子图像
    
    Args:
        mol: RDKit分子对象
        rings: 要高亮的环的列表（每个环是原子索引的元组）
        width: 图像宽度（可选）
        height: 图像高度（可选）
    
    Returns:
        bytes: PNG格式的图像数据
        
    Raises:
        Exception: 图像生成失败时
    """
    try:
        if mol is None:
            raise ValueError("分子对象不能为空")
        
        # 使用配置中的默认尺寸
        img_width = width or MOLECULE_IMAGE_CONFIG['width']
        img_height = height or MOLECULE_IMAGE_CONFIG['height']
        
        # 生成2D坐标（如果没有的话）
        if mol.GetNumConformers() == 0:
            rdDepictor.Compute2DCoords(mol)
        
        # 创建绘图对象
        drawer = rdMolDraw2D.MolDraw2DCairo(img_width, img_height)
        
        # 配置绘图选项
        draw_options = drawer.drawOptions()
        draw_options.addAtomIndices = True  # 显示原子编号
        draw_options.bondLineWidth = MOLECULE_IMAGE_CONFIG['bond_width']
        draw_options.atomLabelFontSize = MOLECULE_IMAGE_CONFIG['atom_label_font_size']
        
        # 设置原子颜色
        for atom_symbol, color in MOLECULE_IMAGE_CONFIG['atom_colors'].items():
            draw_options.updateAtomPalette({atom_symbol: color})
        
        # 准备高亮信息
        atom_highlights = {}
        bond_highlights = {}
        highlight_radii = {}
        
        # 为每个环设置不同的高亮颜色
        colors = MOLECULE_IMAGE_CONFIG['colors']
        for i, ring in enumerate(rings):
            color = colors[i % len(colors)]
            
            # 高亮环中的原子
            for atom_idx in ring:
                if atom_idx < mol.GetNumAtoms():  # 确保原子索引有效
                    atom_highlights[atom_idx] = color
                    highlight_radii[atom_idx] = 0.5
            
            # 高亮环中的键
            for j in range(len(ring)):
                atom1 = ring[j]
                atom2 = ring[(j + 1) % len(ring)]
                
                if atom1 < mol.GetNumAtoms() and atom2 < mol.GetNumAtoms():
                    bond = mol.GetBondBetweenAtoms(atom1, atom2)
                    if bond is not None:
                        bond_highlights[bond.GetIdx()] = color
        
        logger.info(f"高亮 {len(atom_highlights)} 个原子和 {len(bond_highlights)} 个键")
        
        # 绘制分子
        drawer.DrawMolecule(
            mol,
            highlightAtoms=list(atom_highlights.keys()) if atom_highlights else None,
            highlightAtomColors=atom_highlights if atom_highlights else None,
            highlightBonds=list(bond_highlights.keys()) if bond_highlights else None,
            highlightBondColors=bond_highlights if bond_highlights else None,
            highlightAtomRadii=highlight_radii if highlight_radii else None
        )
        
        drawer.FinishDrawing()
        
        # 获取图像数据
        png_data = drawer.GetDrawingText()
        
        if not png_data:
            raise Exception("图像数据为空")
        
        logger.info(f"成功生成分子图像，大小: {len(png_data)} 字节")
        return png_data
        
    except Exception as e:
        logger.error(f"生成分子图像失败: {str(e)}")
        raise Exception(f"图像生成失败: {str(e)}")

def generate_simple_molecule_image(mol: Chem.Mol, 
                                 width: Optional[int] = None, 
                                 height: Optional[int] = None) -> bytes:
    """
    生成简单的分子图像（无高亮）
    
    Args:
        mol: RDKit分子对象
        width: 图像宽度（可选）
        height: 图像高度（可选）
    
    Returns:
        bytes: PNG格式的图像数据
    """
    try:
        if mol is None:
            raise ValueError("分子对象不能为空")
        
        # 使用配置中的默认尺寸
        img_width = width or MOLECULE_IMAGE_CONFIG['width']
        img_height = height or MOLECULE_IMAGE_CONFIG['height']
        
        # 生成2D坐标
        if mol.GetNumConformers() == 0:
            rdDepictor.Compute2DCoords(mol)
        
        # 使用RDKit的简单绘图功能
        img = Draw.MolToImage(mol, size=(img_width, img_height))
        
        # 转换为字节
        import io
        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        png_data = img_bytes.getvalue()
        
        logger.info(f"成功生成简单分子图像，大小: {len(png_data)} 字节")
        return png_data
        
    except Exception as e:
        logger.error(f"生成简单分子图像失败: {str(e)}")
        raise Exception(f"简单图像生成失败: {str(e)}")

def generate_molecule_grid_image(mols: List[Chem.Mol], 
                               legends: Optional[List[str]] = None,
                               mols_per_row: int = 4,
                               sub_img_size: Tuple[int, int] = (200, 200)) -> bytes:
    """
    生成分子网格图像
    
    Args:
        mols: RDKit分子对象列表
        legends: 分子标签列表（可选）
        mols_per_row: 每行分子数量
        sub_img_size: 单个分子图像尺寸
    
    Returns:
        bytes: PNG格式的图像数据
    """
    try:
        if not mols:
            raise ValueError("分子列表不能为空")
        
        # 过滤无效分子
        valid_mols = [mol for mol in mols if mol is not None]
        if not valid_mols:
            raise ValueError("没有有效的分子对象")
        
        # 生成网格图像
        img = Draw.MolsToGridImage(
            valid_mols,
            molsPerRow=mols_per_row,
            subImgSize=sub_img_size,
            legends=legends[:len(valid_mols)] if legends else None
        )
        
        # 转换为字节
        import io
        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        png_data = img_bytes.getvalue()
        
        logger.info(f"成功生成分子网格图像，包含 {len(valid_mols)} 个分子")
        return png_data
        
    except Exception as e:
        logger.error(f"生成分子网格图像失败: {str(e)}")
        raise Exception(f"网格图像生成失败: {str(e)}")

def save_molecule_image(mol: Chem.Mol, filepath: str, 
                       rings: Optional[List[Tuple[int, ...]]] = None,
                       width: Optional[int] = None, 
                       height: Optional[int] = None) -> bool:
    """
    保存分子图像到文件
    
    Args:
        mol: RDKit分子对象
        filepath: 保存路径
        rings: 要高亮的环（可选）
        width: 图像宽度（可选）
        height: 图像高度（可选）
    
    Returns:
        bool: 是否保存成功
    """
    try:
        if rings:
            png_data = generate_molecule_image(mol, rings, width, height)
        else:
            png_data = generate_simple_molecule_image(mol, width, height)
        
        with open(filepath, 'wb') as f:
            f.write(png_data)
        
        logger.info(f"分子图像已保存到: {filepath}")
        return True
        
    except Exception as e:
        logger.error(f"保存分子图像失败: {str(e)}")
        return False

def get_image_config() -> dict:
    """
    获取图像配置信息
    
    Returns:
        dict: 图像配置字典
    """
    return MOLECULE_IMAGE_CONFIG.copy()

def validate_image_parameters(width: Optional[int] = None, 
                            height: Optional[int] = None) -> Tuple[bool, str]:
    """
    验证图像参数
    
    Args:
        width: 图像宽度
        height: 图像高度
    
    Returns:
        Tuple[bool, str]: (是否有效, 错误信息)
    """
    try:
        if width is not None:
            if not isinstance(width, int) or width <= 0:
                return False, "图像宽度必须是正整数"
            if width > 2000:
                return False, "图像宽度不能超过2000像素"
        
        if height is not None:
            if not isinstance(height, int) or height <= 0:
                return False, "图像高度必须是正整数"
            if height > 2000:
                return False, "图像高度不能超过2000像素"
        
        return True, ""
        
    except Exception as e:
        return False, f"参数验证失败: {str(e)}"