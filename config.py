import logging
import os
from functools import lru_cache

# ================ 路径配置 ================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.path.join(BASE_DIR, 'best_gb_model.pkl')
SCALER_PATH = os.path.join(BASE_DIR, 'scaler.pkl')
CACHE_DIR = os.path.join(BASE_DIR, 'cache')
TEMP_DIR = os.path.join(BASE_DIR, 'temp')

# 确保缓存和临时目录存在
os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(TEMP_DIR, exist_ok=True)

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
    'S': 2.58
}

# ================ 性能优化配置 ================
# 缓存大小
CACHE_SIZE = 128  # LRU缓存大小
PARALLEL_PROCESSES = 4  # 并行处理数量

# ================ XTB优化配置 ================
XTB_CONFIG = {
    'opt_level': 'tight',
    'accuracy': 0.1,
    'max_cycles': 1000,
    'max_iterations': 500,
    'electronic_temp': 300,
    'enabled': True,  # 是否启用XTB优化
    'cleanup_temp_files': True  # 是否清理临时文件
}

# ================ 日志配置 ================
LOG_CONFIG = {
    'level': logging.INFO,
    'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    'datefmt': '%Y-%m-%d %H:%M:%S',
    'filename': os.path.join(BASE_DIR, 'app.log'),
    'filemode': 'a'
}

# ================ UI配置 ================
# 分子编辑器配置
JSME_OPTIONS = {
    'width': '600px',
    'height': '400px',
    'options': 'oldlook,star'
}

# 分子图像生成配置
MOLECULE_IMAGE_CONFIG = {
    'width': 400,
    'height': 400,
    'colors': [
        (0.8, 0.1, 0.1),  # 红色
        (0.1, 0.8, 0.1),  # 绿色
        (0.1, 0.1, 0.8),  # 蓝色
        (0.8, 0.8, 0.1),  # 黄色
        (0.8, 0.1, 0.8),  # 紫色
        (0.1, 0.8, 0.8)   # 青色
    ],
    'atom_colors': {
        'C': (0.2, 0.2, 0.2),
        'N': (0.0, 0.0, 1.0),
        'O': (1.0, 0.0, 0.0),
        'S': (1.0, 0.8, 0.0),
        'F': (0.0, 1.0, 0.0),
        'H': (0.8, 0.8, 0.8)
    },
    'highlight_opacity': 0.3,
    'bond_width': 2.0
}

# ================ 缓存装饰器 ================
# 为常用函数提供缓存功能
def cached(func):
    """为函数提供LRU缓存功能的装饰器"""
    return lru_cache(maxsize=CACHE_SIZE)(func)