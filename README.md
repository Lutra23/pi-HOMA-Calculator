# π-HOMA计算器

一个基于Web的分子芳香性计算工具,用于计算分子中芳香环系的π-HOMA (π-Harmonic Oscillator Model of Aromaticity) 值。

## 功能特点

- 🎨 交互式分子结构编辑器
- 🔍 支持SMILES格式输入
- 📊 自动识别和计算芳香环系
- 🖼️ 分子结构可视化
- 🌙 支持深色模式
- 📱 响应式设计,支持移动设备

## 技术栈

- 后端:
  - Python 3.8+
  - Flask - Web框架
  - RDKit - 分子操作和分析
  - scikit-learn - 机器学习模型
  - joblib - 模型序列化
  - XTB (可选,用于分子构型优化)

- 前端:
  - JSME分子编辑器
  - Bootstrap 5
  - Font Awesome
  - JavaScript (原生)

## 安装说明

1. 创建并激活虚拟环境:
```bash
python -m venv venv
source venv/bin/activate  # Linux/Mac
# 或
venv\Scripts\activate  # Windows
```

2. 安装Python依赖:
```bash
pip install -r requirements.txt
```

3. 安装JSME:
```bash
npm install
```

4. 安装XTB (可选):
请参考[XTB官方安装指南](https://xtb-docs.readthedocs.io/en/latest/setup.html)

## 使用说明

1. 启动应用:
```bash
python app.py
```

2. 访问 http://localhost:5000

3. 使用方式:
   - 使用分子编辑器绘制结构
   - 输入SMILES字符串
   - 从预设模板中选择

4. 计算结果包含:
   - 分子2D结构图
   - 每个芳香环的π-HOMA值
   - 环系原子编号

## π-HOMA值说明

π-HOMA (π-Harmonic Oscillator Model of Aromaticity) 是一种评估分子芳香性的指标:

- 取值范围: 0 到 1
- 1 表示完全芳香性
- 0 表示完全非芳香性
- 通常认为 > 0.5 具有显著芳香性

## 开发说明

### 项目结构
```
pi-homa-calculator/
├── app.py              # Flask应用入口
├── config.py           # 配置文件
├── models/             # 核心计算模块
│   ├── calculator.py   # π-HOMA计算器
│   ├── features.py     # 特征提取
│   └── molecule.py     # 分子数据结构
├── utils/              # 工具函数
│   ├── validation.py   # 输入验证
│   └── image.py        # 图像生成
├── static/            # 静态资源
│   ├── css/
│   ├── js/
│   └── jsme/         # JSME编辑器
└── templates/         # HTML模板
```

### 开发环境配置

1. 安装开发依赖:
```bash
pip install -r requirements-dev.txt
```

2. 配置pre-commit hooks:
```bash
pre-commit install
```

3. 运行测试:
```bash
pytest
```

### 开发工具

- **black**: Python代码格式化
- **flake8**: 代码风格检查
- **mypy**: 类型检查
- **pytest**: 单元测试
- **pre-commit**: Git提交前的代码质量检查

## 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

## 致谢

- [JSME Molecule Editor](https://jsme-editor.github.io/)
- [RDKit](https://www.rdkit.org/)
- [XTB](https://xtb-docs.readthedocs.io/)