from flask import Flask, render_template, request, jsonify, send_from_directory
from flask_cors import CORS
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_caching import Cache
import base64
import os
import logging
import time
from datetime import datetime
from models.molecule import MoleculeData
from models.calculator import PiHOMACalculator
from utils.validation import validate_molecule_input, is_aromatic_ring
from utils.image import generate_molecule_image
from config import LOG_CONFIG, CACHE_DIR

# 配置日志
logging.basicConfig(**LOG_CONFIG)
logger = logging.getLogger(__name__)

# 创建Flask应用
app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'dev-secret-key-change-in-production')

# 配置CORS
CORS(app, origins=['http://localhost:5000', 'http://127.0.0.1:5000'])

# 配置缓存
cache_config = {
    'CACHE_TYPE': 'simple',
    'CACHE_DEFAULT_TIMEOUT': 300
}
cache = Cache(app, config=cache_config)

# 配置速率限制
limiter = Limiter(
    app,
    key_func=get_remote_address,
    default_limits=["200 per day", "50 per hour"]
)

# 确保必要的目录存在
os.makedirs(CACHE_DIR, exist_ok=True)

# 添加静态文件路径
@app.route('/static/jsme/<path:filename>')
def serve_jsme(filename):
    """提供JSME文件服务"""
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
    """主页"""
    return render_template('index.html')

@app.route('/health')
def health_check():
    """健康检查端点"""
    return jsonify({
        'status': 'healthy',
        'timestamp': datetime.utcnow().isoformat(),
        'version': '1.0.0'
    })

@app.route('/api/info')
def api_info():
    """API信息端点"""
    return jsonify({
        'name': 'π-HOMA Calculator API',
        'version': '1.0.0',
        'description': 'Calculate π-HOMA values for aromatic molecules',
        'endpoints': {
            'calculate': '/calculate',
            'health': '/health',
            'info': '/api/info'
        }
    })

@app.route('/calculate', methods=['POST'])
@limiter.limit("10 per minute")
@cache.memoize(timeout=300)
def calculate():
    """计算π-HOMA值"""
    start_time = time.time()
    
    try:
        # 获取输入 - 支持JSON和表单两种格式
        if request.is_json:
            input_str = request.json.get('smiles', '')
        else:
            input_str = request.form.get('smiles', '')
        
        if not input_str:
            return jsonify({
                'success': False,
                'error': '请输入SMILES字符串',
                'error_code': 'MISSING_INPUT'
            }), 400
            
        logger.info(f"开始计算π-HOMA，输入: {input_str}")
        
        # 验证输入
        try:
            mol, smiles, input_format = validate_molecule_input(input_str)
        except ValueError as e:
            logger.warning(f"输入验证失败: {str(e)}")
            return jsonify({
                'success': False,
                'error': str(e),
                'error_code': 'INVALID_INPUT'
            }), 400
        
        # 计算π-HOMA
        try:
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
                    'error': '未找到芳香环系',
                    'error_code': 'NO_AROMATIC_RINGS'
                }), 400
            
            # 只计算芳香环的π-HOMA值
            results = calculator.calculate_all_rings_pi_homa(molecule.mol)
            logger.info(f"计算结果: {results}")
            
            # 过滤结果，只保留芳香环的结果
            filtered_results = {}
            for i, idx in enumerate(aromatic_ring_indices):
                if idx in results:
                    filtered_results[i] = results[idx]
            
            logger.info(f"过滤后的结果: {filtered_results}")
            
            if not filtered_results:
                return jsonify({
                    'success': False,
                    'error': '无法计算芳香环的π-HOMA值',
                    'error_code': 'CALCULATION_FAILED'
                }), 500
            
            # 生成带有环高亮的分子图像
            try:
                png_data = generate_molecule_image(molecule.mol, aromatic_rings)
                img_str = base64.b64encode(png_data).decode()
            except Exception as e:
                logger.warning(f"图像生成失败: {str(e)}")
                img_str = ""
            
            # 格式化结果
            rings_data = []
            for ring_idx, pi_homa in filtered_results.items():
                ring_atoms = list(aromatic_rings[ring_idx])
                rings_data.append({
                    'ring_number': ring_idx + 1,
                    'atoms': ring_atoms,
                    'pi_homa': round(pi_homa, 4),
                    'aromaticity_level': get_aromaticity_level(pi_homa)
                })
            
            # 计算总体统计
            pi_homa_values = [ring['pi_homa'] for ring in rings_data]
            statistics = {
                'total_rings': len(rings_data),
                'average_pi_homa': round(sum(pi_homa_values) / len(pi_homa_values), 4),
                'min_pi_homa': min(pi_homa_values),
                'max_pi_homa': max(pi_homa_values)
            }
            
            calculation_time = time.time() - start_time
            
            response_data = {
                'success': True,
                'molecule_image': img_str,
                'rings': rings_data,
                'statistics': statistics,
                'input_format': input_format,
                'smiles': smiles,
                'calculation_time': round(calculation_time, 3),
                'timestamp': datetime.utcnow().isoformat()
            }
            
            logger.info(f"计算完成，耗时: {calculation_time:.3f}秒")
            return jsonify(response_data)
            
        except Exception as e:
            logger.error(f"计算过程出错: {str(e)}", exc_info=True)
            return jsonify({
                'success': False,
                'error': '计算过程中发生错误',
                'error_code': 'CALCULATION_ERROR',
                'details': str(e) if app.debug else None
            }), 500
        
    except Exception as e:
        logger.error(f"请求处理出错: {str(e)}", exc_info=True)
        return jsonify({
            'success': False,
            'error': '服务器内部错误',
            'error_code': 'INTERNAL_ERROR',
            'details': str(e) if app.debug else None
        }), 500

def get_aromaticity_level(pi_homa):
    """根据π-HOMA值判断芳香性水平"""
    if pi_homa >= 0.8:
        return 'strong'
    elif pi_homa >= 0.5:
        return 'moderate'
    elif pi_homa >= 0.2:
        return 'weak'
    else:
        return 'non-aromatic'

@app.errorhandler(404)
def not_found(error):
    """404错误处理"""
    return jsonify({
        'success': False,
        'error': '页面未找到',
        'error_code': 'NOT_FOUND'
    }), 404

@app.errorhandler(500)
def internal_error(error):
    """500错误处理"""
    logger.error(f"内部服务器错误: {str(error)}")
    return jsonify({
        'success': False,
        'error': '服务器内部错误',
        'error_code': 'INTERNAL_ERROR'
    }), 500

@app.errorhandler(429)
def ratelimit_handler(e):
    """速率限制错误处理"""
    return jsonify({
        'success': False,
        'error': '请求过于频繁，请稍后重试',
        'error_code': 'RATE_LIMIT_EXCEEDED'
    }), 429

# 添加安全头
@app.after_request
def after_request(response):
    """添加安全头"""
    response.headers['X-Content-Type-Options'] = 'nosniff'
    response.headers['X-Frame-Options'] = 'DENY'
    response.headers['X-XSS-Protection'] = '1; mode=block'
    response.headers['Strict-Transport-Security'] = 'max-age=31536000; includeSubDomains'
    return response

if __name__ == '__main__':
    # 开发环境配置
    debug_mode = os.environ.get('FLASK_DEBUG', 'False').lower() == 'true'
    port = int(os.environ.get('PORT', 5000))
    host = os.environ.get('HOST', '0.0.0.0')
    
    logger.info(f"启动π-HOMA计算器服务器 - Debug: {debug_mode}, Host: {host}, Port: {port}")
    
    app.run(
        host=host,
        port=port,
        debug=debug_mode,
        threaded=True
    )