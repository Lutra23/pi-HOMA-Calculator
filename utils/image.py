"""分子图像生成工具"""
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor
from config import MOLECULE_IMAGE_CONFIG

def generate_molecule_image(mol, rings):
    """
    生成带有环高亮的分子图像
    
    Args:
        mol: RDKit分子对象
        rings: 要高亮的环的列表（每个环是原子索引的列表）
    
    Returns:
        bytes: PNG格式的图像数据
    """
    # 生成2D坐标（如果没有的话）
    rdDepictor.Compute2DCoords(mol)
    
    # 创建绘图对象
    drawer = Draw.rdMolDraw2D.MolDraw2DCairo(
        MOLECULE_IMAGE_CONFIG['width'],
        MOLECULE_IMAGE_CONFIG['height']
    )
    drawer.drawOptions().addAtomIndices = True  # 显示原子编号
    
    # 准备高亮信息
    atom_highlights = {}
    bond_highlights = {}
    highlight_radii = {}
    
    # 为每个环设置高亮
    for i, ring in enumerate(rings):
        color = MOLECULE_IMAGE_CONFIG['colors'][i % len(MOLECULE_IMAGE_CONFIG['colors'])]
        # 高亮原子
        for atom_idx in ring:
            atom_highlights[atom_idx] = color
            highlight_radii[atom_idx] = 0.5
        # 高亮键
        for j in range(len(ring)):
            atom1 = ring[j]
            atom2 = ring[(j + 1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(atom1, atom2)
            if bond is not None:
                bond_highlights[bond.GetIdx()] = color
    
    # 绘制分子
    drawer.DrawMolecule(
        mol,
        highlightAtoms=list(atom_highlights.keys()),
        highlightAtomColors=atom_highlights,
        highlightBonds=list(bond_highlights.keys()),
        highlightBondColors=bond_highlights,
        highlightAtomRadii=highlight_radii
    )
    drawer.FinishDrawing()
    
    # 获取图像数据
    png_data = drawer.GetDrawingText()
    return png_data