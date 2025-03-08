from flask import Flask, render_template, request, jsonify, send_from_directory
import base64
from models.molecule import MoleculeData
from models.calculator import PiHOMACalculator
from utils.validation import validate_molecule_input, is_aromatic_ring
from utils.image import generate_molecule_image
import logging
import os
from config import LOG_CONFIG

# 配置日志
logging.basicConfig(**LOG_CONFIG)
logger = logging.getLogger(__name__)

app = Flask(__name__)

# 添加静态文件路径
@app.route('/static/jsme/<path:filename>')
def serve_jsme(filename):
    logger.info(f'请求JSME文件: {filename}')
    jsme_dir = os.path.join(app.static_folder, 'jsme')
    if os.path.exists(os.path.join(jsme_dir, filename)):
        logger.info(f'找到文件: {filename}')
        return send_from_directory(jsme_dir, filename)
    else:
        logger.error(f'文件未找到: {filename}')
        return f'File not found: {filename}', 404

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/calculate', methods=['POST'])
def calculate():
    try:
        # 获取输入 - 支持JSON和表单两种格式
        if request.is_json:
            input_str = request.json.get('smiles', '')
        else:
            input_str = request.form.get('smiles', '')
        
        if not input_str:
            return jsonify({
                'success': False,
                'error': '请输入SMILES字符串'
            })
            
        # 验证输入
        mol, smiles, input_format = validate_molecule_input(input_str)
        
        # 计算π-HOMA
        molecule = MoleculeData(smiles)
        calculator = PiHOMACalculator()
        
        # 获取所有环
        all_rings = list(molecule.mol.GetRingInfo().AtomRings())
        logger.info(f"找到的环: {all_rings}")
        
        # 筛选芳香环
        aromatic_rings = []
        aromatic_ring_indices = []
        for i, ring in enumerate(all_rings):
            if is_aromatic_ring(molecule.mol, ring):
                aromatic_rings.append(ring)
                aromatic_ring_indices.append(i)
        
        logger.info(f"芳香环: {aromatic_rings}")
        
        if not aromatic_rings:
            return jsonify({
                'success': False,
                'error': '未找到芳香环系'
            })
        
        # 只计算芳香环的π-HOMA值
        results = calculator.calculate_all_rings_pi_homa(molecule.mol)
        logger.info(f"计算结果: {results}")
        filtered_results = {i: results[idx] for i, idx in enumerate(aromatic_ring_indices)}
        logger.info(f"过滤后的结果: {filtered_results}")
        
        # 生成带有环高亮的分子图像
        png_data = generate_molecule_image(molecule.mol, aromatic_rings)
        img_str = base64.b64encode(png_data).decode()
        
        # 格式化结果
        rings_data = []
        for ring_idx, pi_homa in filtered_results.items():
            ring_atoms = list(aromatic_rings[ring_idx])
            rings_data.append({
                'ring_number': ring_idx + 1,
                'atoms': ring_atoms,
                'pi_homa': round(pi_homa, 4),
                'format': input_format
            })
        
        response_data = {
            'success': True,
            'molecule_image': img_str,
            'rings': rings_data
        }
        logger.info(f"返回数据结构: {response_data}")
        
        return jsonify(response_data)
        
    except Exception as e:
        logger.error(f"计算过程出错: {str(e)}")
        return jsonify({
            'success': False,
            'error': str(e)
        })

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)