# 🚀 π-HOMA-Calculator: A Web-Based Tool for Calculating Molecular Aromaticity 🌟

## About the Project

The π-HOMA-Calculator is a web-based tool for calculating the π-HOMA (π-Harmonic Oscillator Model of Aromaticity) values of molecules. This tool is designed to provide a user-friendly interface for calculating and visualizing π-HOMA values.

## Features

### 🎨 Interactive Molecular Structure Editor

The π-HOMA-Calculator comes with an interactive molecular structure editor, allowing users to draw molecular structures and calculate π-HOMA values.

### 🔍 SMILES Format Input

Users can input molecular structures in SMILES format, allowing for easy import and calculation of π-HOMA values.

### 📊 Automatic Identification and Calculation of Aromatic Rings

The π-HOMA-Calculator automatically identifies and calculates π-HOMA values for aromatic rings in the input molecular structure.

### 🖼️ Molecular Structure Visualization

The π-HOMA-Calculator provides interactive molecular structure visualization, allowing users to examine the molecular structure and π-HOMA values.

### 🌙 Dark Mode Support

The π-HOMA-Calculator supports dark mode, providing a visually appealing interface for users who prefer a darker theme.

### 📱 Responsive Design for Mobile Devices

The π-HOMA-Calculator is optimized for mobile devices, providing a responsive design that adapts to different screen sizes and devices.

## Technical Details

### 👥 Frontend

- JavaScript (native)
- Bootstrap 5
- Font Awesome
- JSME Molecule Editor

### 🕳️ Backend

- Python 3.8+
- Flask - Web framework
- RDKit - molecular operation and analysis
- scikit-learn - machine learning model
- joblib - model serialization
- XTB (optional) - molecular conformation optimization

## Installation and Setup

### 1️⃣ Create and Activate a Virtual Environment

```bash
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or
venv\Scripts\activate  # Windows
```

### 2️⃣ Install Python Dependencies

```bash
pip install -r requirements.txt
```

### 3️⃣ Install JSME

```bash
npm install
```

### 4️⃣ Install XTB (optional)

Please refer to the [XTB official installation guide](https://xtb-docs.readthedocs.io/en/latest/setup.html)

## Usage

### 1️⃣ Run the Application

```bash
python app.py
```

### 2️⃣ Access the π-HOMA-Calculator

[http://localhost:5000](http://localhost:5000)

### 3️⃣ Use the π-HOMA-Calculator

- Use the molecular editor to draw molecular structures
- Input SMILES strings
- Select from pre-defined templates

### 4️⃣ View Results

The π-HOMA-Calculator provides the following results:
- 2D molecular structure diagram
- π-HOMA values for each aromatic ring
- Ring atom numbering

## π-HOMA Value Explanation

π-HOMA values are a measure of a molecule's aromaticity:

- Values range from 0 to 1
- 0 indicates no aromaticity
- 1 indicates complete aromaticity
- Values > 0.5 are typically considered to exhibit significant aromaticity

## Contributing Guidelines

### 📸 Project Structure

The π-HOMA-Calculator repository is organized as follows:
```
pi-homa-calculator/
├── app.py              # Flask application entry point
├── config.py           # configuration file
├── models/             # core calculation modules
│   ├── calculator.py   # π-HOMA calculator
│   ├── features.py     # feature extraction
│   └── molecule.py     # molecular data structure
├── utils/              # utility functions
│   ├── validation.py   # input validation
│   └── image.py        # image generation
├── static/            # static resources
│   ├── css/
│   ├── js/
│   └──yme/         # JSME editor
└── templates/         # HTML templates
```

### 🔧 Development Environment Setup

1. Install development dependencies:
```bash
pip install -r requirements-dev.txt
```

2. Configure pre-commit hooks:
```bash
pre-commit install
```

3. Run tests:
```bash
pytest
```

### 💡 Development Tools

- **black**: Python code formatter
- **flake8**: code style checker
- **mypy**: type checker
- **pytest**: unit testing
- **pre-commit**: Git commit quality checker

## 😊 License

The π-HOMA-Calculator is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- [JSME Molecule Editor](https://jsme-editor.github.io/)
- [RDKit](https://www.rdkit.org/)
- [XTB](https://xtb-docs.readthedocs.io/)
