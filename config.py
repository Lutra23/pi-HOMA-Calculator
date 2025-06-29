import logging
import os
from functools import lru_cache
from pathlib import Path

# ================ 路径配置 ================
BASE_DIR = Path(__file__).parent.absolute()
MODEL_PATH = BASE_DIR / 'best_gb_model.pkl'
SCALER_PATH = BASE_DIR / 'scaler.pkl'
CACHE_DIR = BASE_DIR / 'cache'
TEMP_DIR = BASE_DIR / 'temp'
LOGS_DIR = BASE_DIR / 'logs'

# 确保必要目录存在
for directory in [CACHE_DIR, TEMP_DIR, LOGS_DIR]:
    directory.mkdir(exist_ok=True)

# ================ 环境配置 ================
ENVIRONMENT = os.getenv('ENVIRONMENT', 'development')
DEBUG = os.getenv('DEBUG', 'False').lower() == 'true'
SECRET_KEY = os.getenv('SECRET_KEY', 'dev-secret-key-change-in-production')

# ================ 分子计算配置 ================
# 理想键级值配置
IDEAL_BOND_ORDERS = {
    ('C', 'N'): 0.5129,
    ('N', 'C'): 0.5129,
    ('C', 'S'): 0.4497,
    ('S', 'C'): 0.4497,
    ('C', 'O'): 0.43615,
    ('O', 'C'): 0.43615,
    ('C', 'C'): 0.66,
}

# 电负性配置
ELECTRONEGATIVITY = {
    'H': 2.20,
    'C': 2.55,
    'N': 3.04,
    'O': 3.44,
    'F': 3.98,
    'S': 2.58,
    'P': 2.19,
    'Cl': 3.16,
    'Br': 2.96,
    'I': 2.66
}

# ================ 性能优化配置 ================
# 缓存大小
CACHE_SIZE = 256  # LRU缓存大小
PARALLEL_PROCESSES = min(4, os.cpu_count() or 1)  # 并行处理数量

# 缓存过期时间（秒）
CACHE_TIMEOUT = {
    'calculation': 3600,  # 计算结果缓存1小时
    'molecule': 1800,     # 分子数据缓存30分钟
    'image': 7200         # 图像缓存2小时
}

# ================ XTB优化配置 ================
XTB_CONFIG = {
    'opt_level': 'tight',
    'accuracy': 0.1,
    'max_cycles': 1000,
    'max_iterations': 500,
    'electronic_temp': 300,
    'enabled': True,  # 是否启用XTB优化
    'cleanup_temp_files': True,  # 是否清理临时文件
    'timeout': 300  # XTB计算超时时间（秒）
}

# ================ 日志配置 ================
LOG_LEVEL = getattr(logging, os.getenv('LOG_LEVEL', 'INFO').upper())

LOG_CONFIG = {
    'level': LOG_LEVEL,
    'format': '%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s',
    'datefmt': '%Y-%m-%d %H:%M:%S',
    'handlers': [
        logging.StreamHandler(),
        logging.FileHandler(LOGS_DIR / 'app.log', encoding='utf-8')
    ]
}

# 配置日志轮转
if ENVIRONMENT == 'production':
    from logging.handlers import RotatingFileHandler
    LOG_CONFIG['handlers'] = [
        logging.StreamHandler(),
        RotatingFileHandler(
            LOGS_DIR / 'app.log',
            maxBytes=10*1024*1024,  # 10MB
            backupCount=5,
            encoding='utf-8'
        )
    ]

# ================ Web服务配置 ================
# Flask配置
FLASK_CONFIG = {
    'SECRET_KEY': SECRET_KEY,
    'JSON_AS_ASCII': False,
    'JSONIFY_PRETTYPRINT_REGULAR': DEBUG,
    'MAX_CONTENT_LENGTH': 16 * 1024 * 1024,  # 16MB
}

# CORS配置
CORS_CONFIG = {
    'origins': os.getenv('CORS_ORIGINS', 'http://localhost:*,http://127.0.0.1:*').split(','),
    'methods': ['GET', 'POST', 'OPTIONS'],
    'allow_headers': ['Content-Type', 'Authorization']
}

# 速率限制配置
RATE_LIMIT_CONFIG = {
    'default': "200 per day, 50 per hour",
    'calculate': "30 per minute, 10 per second",
    'storage_uri': os.getenv('REDIS_URL', 'memory://')
}

# ================ UI配置 ================
# 分子编辑器配置
JSME_OPTIONS = {
    'width': '100%',
    'height': '400px',
    'options': 'oldlook,star,polarnitro,nocanonize,nostereo'
}

# 分子图像生成配置
MOLECULE_IMAGE_CONFIG = {
    'width': 500,
    'height': 400,
    'dpi': 150,
    'colors': [
        (0.8, 0.1, 0.1),  # 红色
        (0.1, 0.8, 0.1),  # 绿色
        (0.1, 0.1, 0.8),  # 蓝色
        (0.8, 0.8, 0.1),  # 黄色
        (0.8, 0.1, 0.8),  # 紫色
        (0.1, 0.8, 0.8),  # 青色
        (0.8, 0.5, 0.1),  # 橙色
        (0.5, 0.1, 0.8),  # 紫罗兰
    ],
    'atom_colors': {
        'C': (0.2, 0.2, 0.2),
        'N': (0.0, 0.0, 1.0),
        'O': (1.0, 0.0, 0.0),
        'S': (1.0, 0.8, 0.0),
        'F': (0.0, 1.0, 0.0),
        'H': (0.8, 0.8, 0.8),
        'P': (1.0, 0.5, 0.0),
        'Cl': (0.0, 1.0, 0.0),
        'Br': (0.5, 0.2, 0.1),
        'I': (0.5, 0.0, 0.5)
    },
    'highlight_opacity': 0.3,
    'bond_width': 2.0,
    'atom_label_font_size': 12
}

# ================ 数据库配置（未来扩展用） ================
DATABASE_CONFIG = {
    'url': os.getenv('DATABASE_URL', f'sqlite:///{BASE_DIR}/app.db'),
    'echo': DEBUG,
    'pool_size': 10,
    'max_overflow': 20,
    'pool_timeout': 30,
    'pool_recycle': 3600
}

# ================ 安全配置 ================
SECURITY_CONFIG = {
    'max_smiles_length': 1000,
    'allowed_atoms': {'C', 'N', 'O', 'S', 'F', 'H', 'P', 'Cl', 'Br', 'I'},
    'max_atoms': 200,
    'max_rings': 20,
    'session_timeout': 3600,  # 1小时
    'csrf_protection': ENVIRONMENT == 'production'
}

# ================ 监控配置 ================
MONITORING_CONFIG = {
    'enable_metrics': ENVIRONMENT == 'production',
    'metrics_endpoint': '/metrics',
    'health_check_endpoint': '/health',
    'performance_tracking': True
}

# ================ 缓存装饰器 ================
def cached(func):
    """为函数提供LRU缓存功能的装饰器"""
    return lru_cache(maxsize=CACHE_SIZE)(func)

# ================ 验证配置 ================
def validate_config():
    """验证配置的有效性"""
    errors = []
    
    # 检查必要文件
    if not MODEL_PATH.exists():
        errors.append(f"模型文件不存在: {MODEL_PATH}")
    
    if not SCALER_PATH.exists():
        errors.append(f"标准化器文件不存在: {SCALER_PATH}")
    
    # 检查环境变量
    if ENVIRONMENT == 'production' and SECRET_KEY == 'dev-secret-key-change-in-production':
        errors.append("生产环境必须设置安全的SECRET_KEY")
    
    # 检查数值配置
    if CACHE_SIZE <= 0:
        errors.append("CACHE_SIZE必须大于0")
    
    if PARALLEL_PROCESSES <= 0:
        errors.append("PARALLEL_PROCESSES必须大于0")
    
    if errors:
        raise ValueError("配置验证失败:\n" + "\n".join(errors))

# 在导入时验证配置
if __name__ != '__main__':
    try:
        validate_config()
    except ValueError as e:
        print(f"配置错误: {e}")
        if ENVIRONMENT == 'production':
            raise

# ================ 导出配置 ================
__all__ = [
    'BASE_DIR', 'MODEL_PATH', 'SCALER_PATH', 'CACHE_DIR', 'TEMP_DIR', 'LOGS_DIR',
    'ENVIRONMENT', 'DEBUG', 'SECRET_KEY',
    'IDEAL_BOND_ORDERS', 'ELECTRONEGATIVITY',
    'CACHE_SIZE', 'PARALLEL_PROCESSES', 'CACHE_TIMEOUT',
    'XTB_CONFIG', 'LOG_CONFIG',
    'FLASK_CONFIG', 'CORS_CONFIG', 'RATE_LIMIT_CONFIG',
    'JSME_OPTIONS', 'MOLECULE_IMAGE_CONFIG',
    'DATABASE_CONFIG', 'SECURITY_CONFIG', 'MONITORING_CONFIG',
    'cached', 'validate_config'
]