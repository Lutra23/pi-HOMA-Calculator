// 分子编辑器相关功能
class MoleculeEditor {
    constructor() {
        this.jsmeApplet = null;
        this.isInitialized = false;
        this.initEditor();
    }

    initEditor() {
        // 等待JSME加载完成
        if (typeof JSApplet !== 'undefined') {
            this.setupJSME();
        } else {
            // 如果JSME还没加载，等待一段时间后重试
            setTimeout(() => this.initEditor(), 500);
        }
    }

    setupJSME() {
        try {
            const container = document.getElementById('jsme_container');
            if (!container) {
                console.error('JSME容器未找到');
                return;
            }

            // JSME配置
            const jsmeOptions = {
                'options': 'oldlook,star,polarnitro,nocanonize,nostereo',
                'width': '100%',
                'height': '400px'
            };

            // 创建JSME实例
            this.jsmeApplet = new JSApplet.JSME('jsme_container', '600px', '400px', jsmeOptions);
            
            if (this.jsmeApplet) {
                this.isInitialized = true;
                this.setupEventListeners();
                console.log('JSME编辑器初始化成功');
            }
        } catch (error) {
            console.error('JSME初始化失败:', error);
            this.showFallbackEditor();
        }
    }

    setupEventListeners() {
        if (!this.jsmeApplet) return;

        // 监听分子变化事件
        try {
            this.jsmeApplet.setCallBack('AfterStructureModified', () => {
                this.onMoleculeChanged();
            });
        } catch (error) {
            console.warn('无法设置JSME回调:', error);
        }
    }

    onMoleculeChanged() {
        // 分子结构改变时的处理
        if (this.jsmeApplet) {
            const smiles = this.getSMILES();
            if (smiles) {
                // 可以在这里添加实时验证或其他功能
                this.validateMolecule(smiles);
            }
        }
    }

    getSMILES() {
        if (!this.jsmeApplet || !this.isInitialized) {
            return '';
        }

        try {
            return this.jsmeApplet.smiles() || '';
        } catch (error) {
            console.error('获取SMILES失败:', error);
            return '';
        }
    }

    setSMILES(smiles) {
        if (!this.jsmeApplet || !this.isInitialized) {
            console.warn('JSME编辑器未初始化');
            return false;
        }

        try {
            this.jsmeApplet.readGenericMolecularInput(smiles);
            return true;
        } catch (error) {
            console.error('设置SMILES失败:', error);
            return false;
        }
    }

    clear() {
        if (!this.jsmeApplet || !this.isInitialized) {
            return;
        }

        try {
            this.jsmeApplet.clear();
        } catch (error) {
            console.error('清空编辑器失败:', error);
        }
    }

    undo() {
        if (!this.jsmeApplet || !this.isInitialized) {
            return;
        }

        try {
            // JSME可能没有直接的undo方法，这里需要根据实际API调整
            console.log('撤销操作');
        } catch (error) {
            console.error('撤销操作失败:', error);
        }
    }

    redo() {
        if (!this.jsmeApplet || !this.isInitialized) {
            return;
        }

        try {
            // JSME可能没有直接的redo方法，这里需要根据实际API调整
            console.log('重做操作');
        } catch (error) {
            console.error('重做操作失败:', error);
        }
    }

    exportMolecule() {
        const smiles = this.getSMILES();
        if (!smiles) {
            calculator?.showError('编辑器中没有分子结构');
            return;
        }

        // 创建下载链接
        const blob = new Blob([smiles], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = 'molecule.smi';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(url);

        calculator?.showSuccess('分子结构已导出');
    }

    validateMolecule(smiles) {
        // 简单的SMILES验证
        if (!smiles || smiles.trim() === '') {
            return false;
        }

        // 检查是否包含芳香环标记
        const hasAromaticRing = /[a-z]/.test(smiles) || smiles.includes('1') || smiles.includes('2');
        
        if (!hasAromaticRing) {
            // 可以显示警告，但不阻止计算
            console.warn('分子可能不包含芳香环');
        }

        return true;
    }

    showFallbackEditor() {
        const container = document.getElementById('jsme_container');
        if (container) {
            container.innerHTML = `
                <div class="alert alert-warning text-center">
                    <h5><i class="fas fa-exclamation-triangle"></i> 分子编辑器加载失败</h5>
                    <p>请使用SMILES输入方式或常用模板功能</p>
                    <button class="btn btn-primary" onclick="document.getElementById('smiles-tab').click()">
                        切换到SMILES输入
                    </button>
                </div>
            `;
        }
    }

    // 静态方法，用于检查JSME是否可用
    static isJSMEAvailable() {
        return typeof JSApplet !== 'undefined';
    }
}

// 全局编辑器实例
let moleculeEditor;

// 全局函数（保持向后兼容）
function clearEditor() {
    if (moleculeEditor) {
        moleculeEditor.clear();
    }
}

function undoLastAction() {
    if (moleculeEditor) {
        moleculeEditor.undo();
    }
}

function redoLastAction() {
    if (moleculeEditor) {
        moleculeEditor.redo();
    }
}

function exportMolecule() {
    if (moleculeEditor) {
        moleculeEditor.exportMolecule();
    }
}

// 初始化编辑器
document.addEventListener('DOMContentLoaded', () => {
    // 延迟初始化，确保JSME库已加载
    setTimeout(() => {
        moleculeEditor = new MoleculeEditor();
        
        // 将编辑器实例暴露给全局，以便其他模块使用
        window.jsmeApplet = moleculeEditor?.jsmeApplet;
    }, 1000);
});