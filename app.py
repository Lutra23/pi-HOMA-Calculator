from flask import Flask, render_template, request, jsonify, send_from_directory
import base64
from models.molecule import MoleculeData
from models.calculator import PiHOMACalculator
from models.features import MolecularFeatureExtractor # Added import
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
            
        # 验证输入 (This might still raise an error for non-aromatic molecules based on current utils/validation.py)
        # The subtask assumes this validation will be relaxed or handled appropriately for non-aromatic cases.
        mol_obj, smiles, input_format = validate_molecule_input(input_str)
        
        molecule_data = MoleculeData(smiles) # RDKit Mol object is in molecule_data.mol

        # ---- START Molecular Descriptors Calculation ----
        feature_extractor = MolecularFeatureExtractor()
        # Ensure molecule_data.mol is used, which should be the same as mol_obj from validate_molecule_input
        features, feature_names = feature_extractor.get_molecular_descriptors(molecule_data.mol) 
        
        descriptors_map = {
            'mol_weight': 'Molecular Weight',
            'logp': 'LogP',
            'tpsa': 'TPSA',
            'num_rotatable_bonds': 'Number of Rotatable Bonds',
            'num_h_donors': 'Number of H-bond Donors',
            'num_h_acceptors': 'Number of H-bond Acceptors',
            'exact_mol_weight': 'Exact Molecular Weight'
        }
        
        molecular_descriptors_data = {}
        if features and feature_names:
            feature_dict = dict(zip(feature_names, features))
            for internal_name, display_name in descriptors_map.items():
                molecular_descriptors_data[display_name] = feature_dict.get(internal_name, 'N/A') 
        else:
            logger.warning(f"Could not retrieve molecular descriptors for SMILES: {smiles}")
            for display_name in descriptors_map.values():
                 molecular_descriptors_data[display_name] = 'N/A'
        # ---- END Molecular Descriptors Calculation ----

        # Initialize response_data with descriptors
        response_data = {
            'success': True, # Assume success, can be overridden by HOMA errors
            'molecular_descriptors': molecular_descriptors_data,
            'rings': [],
            'molecule_image': None,
            'homa_message': '', # For non-critical HOMA messages
            'homa_error': ''   # For critical HOMA errors
        }

        # ---- START HOMA Calculation Block ----
        try:
            calculator = PiHOMACalculator() # This might raise .pkl error
            
            all_rings = list(molecule_data.mol.GetRingInfo().AtomRings())
            logger.info(f"找到的环: {all_rings}")
            
            aromatic_rings = []
            aromatic_ring_indices = []
            for i, ring in enumerate(all_rings):
                if is_aromatic_ring(molecule_data.mol, ring): # is_aromatic_ring is from utils.validation
                    aromatic_rings.append(ring)
                    aromatic_ring_indices.append(i)
            
            logger.info(f"芳香环: {aromatic_rings}")
            
            if not aromatic_rings:
                response_data['homa_message'] = 'No aromatic rings found for HOMA analysis.'
                # success remains true as descriptors are present
            else:
                homa_results = calculator.calculate_all_rings_pi_homa(molecule_data.mol)
                logger.info(f"HOMA 计算结果: {homa_results}")
                filtered_homa_results = {i: homa_results[idx] for i, idx in enumerate(aromatic_ring_indices)}
                logger.info(f"过滤后的 HOMA 结果: {filtered_homa_results}")
                
                png_data = generate_molecule_image(molecule_data.mol, aromatic_rings)
                img_str = base64.b64encode(png_data).decode()
                response_data['molecule_image'] = img_str
                
                rings_data_list = []
                for ring_idx, pi_homa_value in filtered_homa_results.items():
                    ring_atoms_list = list(aromatic_rings[ring_idx])
                    rings_data_list.append({
                        'ring_number': ring_idx + 1,
                        'atoms': ring_atoms_list,
                        'pi_homa': round(pi_homa_value, 4),
                        'format': input_format 
                    })
                response_data['rings'] = rings_data_list
                response_data['homa_message'] = 'HOMA analysis completed.'

        except Exception as homa_exc:
            logger.error(f"HOMA 计算过程出错: {str(homa_exc)}")
            response_data['homa_error'] = f'HOMA calculation failed: {str(homa_exc)}'
            # Depending on policy, you might set response_data['success'] = False here
            # For now, if descriptors are there, let's keep success = True unless a top-level error occurs
        # ---- END HOMA Calculation Block ----
        
        logger.info(f"最终返回数据结构: {response_data}")
        return jsonify(response_data)
        
    except ValueError as ve: # Catch validation errors from validate_molecule_input specifically
        logger.error(f"输入验证错误: {str(ve)}")
        return jsonify({
            'success': False,
            'error': str(ve),
            'molecular_descriptors': {}, # Ensure structure even on validation error
            'rings': [],
            'molecule_image': None
        })
    except Exception as e:
        logger.error(f"计算过程出现意外错误: {str(e)}")
        return jsonify({
            'success': False,
            'error': f"An unexpected error occurred: {str(e)}",
            'molecular_descriptors': {}, # Ensure structure even on unexpected error
            'rings': [],
            'molecule_image': None
        })

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)