/**
 * 分子编辑器相关的功能实现
 */

// 全局变量存储JSME实例
window.jsmeApplet = null;

/**
 * JSME分子编辑器初始化配置
 */
const JSME_OPTIONS = {
    "options": "nosearch,nopaste,oldlook",
    "borderColor": "#dee2e6",
    "backgroundColor": "#ffffff"
};

// 添加脚本加载状态检查
document.addEventListener('DOMContentLoaded', function() {
    console.log('DOM加载完成,检查JSME状态...');
    
    // 获取JSME脚本元素
    const jsmeScript = document.querySelector('script[src*="jsme.nocache.js"]');
    if (!jsmeScript) {
        console.error('未找到JSME脚本标签');
        showError('JSME脚本未正确加载');
        return;
    }
    
    // 检查JSME是否已加载
    if (typeof JSApplet !== 'undefined') {
        console.log('JSME已加载,直接初始化');
        initializeJSME();
    } else {
        console.log('JSME尚未加载,等待加载完成...');
        // 设置加载超时
        setTimeout(() => {
            if (typeof JSApplet === 'undefined') {
                console.error('JSME加载超时');
                showError('JSME加载超时,请刷新页面重试');
            }
        }, 10000);
    }
});

/**
 * 显示错误信息
 */
function showError(message) {
    console.error(message);
    window.showToast(message, 'error');
    
    const container = document.getElementById('jsme_container');
    if (container) {
        container.innerHTML = `
            <div class="alert alert-danger m-3">
                <h5>分子编辑器加载失败</h5>
                <p>${message}</p>
                <button class="btn btn-outline-danger btn-sm" onclick="location.reload()">
                    <i class="fas fa-sync-alt me-1"></i>刷新页面
                </button>
            </div>
        `;
    }
}

/**
 * 初始化JSME编辑器
 */
function initializeJSME() {
    console.log('开始初始化JSME...');
    
    try {
        // 检查容器元素
        const container = document.getElementById('jsme_container');
        if (!container) {
            throw new Error('找不到JSME容器元素(#jsme_container)');
        }
        console.log('找到JSME容器元素');

        // 检查JSME是否正确加载
        if (typeof JSApplet === 'undefined' || typeof JSApplet.JSME === 'undefined') {
            throw new Error('JSME库未正确加载');
        }
        console.log('JSME库加载正常');

        // 创建编辑器实例
        window.jsmeApplet = new JSApplet.JSME("jsme_container", "100%", "400px", {
            ...JSME_OPTIONS
        });
        console.log('JSME实例创建成功');

        // 初始化分子模板
        initializeMoleculeTemplates();
        console.log('分子模板初始化完成');
        
        // 设置深色模式监听
        observeDarkModeChanges();
        console.log('深色模式监听设置完成');
        
        // 初始化成功后设置主题
        updateEditorTheme();
        console.log('主题设置完成');

        // 通知用户
        window.showToast('分子编辑器加载完成', 'success');
        
    } catch (error) {
        console.error('分子编辑器初始化失败:', error);
        showError(error.message);
    }
}

/**
 * JSME加载完成后的回调函数
 */
window.jsmeOnLoad = function() {
    console.log('JSME加载完成回调触发');
    initializeJSME();
}

// 常用分子模板定义
const MOLECULE_TEMPLATES = [
    {
        name: '苯',
        smiles: 'c1ccccc1',
        description: '单环芳香烃'
    },
    {
        name: '吡啶',
        smiles: 'n1ccccc1',
        description: '含氮杂环'
    },
    {
        name: '噻吩',
        smiles: 's1cccc1',
        description: '含硫杂环'
    },
    {
        name: '呋喃',
        smiles: 'o1cccc1',
        description: '含氧杂环'
    },
    {
        name: '萘',
        smiles: 'c1ccc2ccccc2c1',
        description: '双环芳香烃'
    },
    {
        name: '吲哚',
        smiles: 'c1ccc2c(c1)cc[nH]2',
        description: '含氮稠环'
    }
];

/**
 * 显示提示消息
 * @param {string} message - 消息内容
 * @param {string} type - 消息类型（success/warning/error）
 */
window.showToast = function(message, type = 'info') {
    // 创建toast容器（如果不存在）
    let container = document.querySelector('.toast-container');
    if (!container) {
        container = document.createElement('div');
        container.className = 'toast-container position-fixed bottom-0 end-0 p-3';
        document.body.appendChild(container);
    }

    // 创建新的toast元素
    const toast = document.createElement('div');
    toast.className = `toast align-items-center text-white bg-${type} border-0`;
    toast.setAttribute('role', 'alert');
    toast.setAttribute('aria-live', 'assertive');
    toast.setAttribute('aria-atomic', 'true');
    
    toast.innerHTML = `
        <div class="d-flex">
            <div class="toast-body">
                ${message}
            </div>
            <button type="button" class="btn-close btn-close-white me-2 m-auto" data-bs-dismiss="toast"></button>
        </div>
    `;
    
    // 添加到容器
    container.appendChild(toast);
    
    // 显示toast
    const bsToast = new bootstrap.Toast(toast, {
        autohide: true,
        delay: 3000
    });
    bsToast.show();
    
    // 移除已关闭的toast
    toast.addEventListener('hidden.bs.toast', () => {
        toast.remove();
        if (container.children.length === 0) {
            container.remove();
        }
    });
}

/**
 * 初始化分子模板列表
 */
function initializeMoleculeTemplates() {
    const templatesGrid = document.querySelector('.templates-grid');
    if (!templatesGrid) return;

    MOLECULE_TEMPLATES.forEach(template => {
        const templateCard = document.createElement('div');
        templateCard.className = 'template-card';
        templateCard.innerHTML = `
            <h6>${template.name}</h6>
            <p class="text-muted small mb-2">${template.description}</p>
            <small class="text-secondary">${template.smiles}</small>
        `;
        
        templateCard.addEventListener('click', () => {
            loadMoleculeTemplate(template.smiles);
            // 切换到编辑器标签
            document.getElementById('editor-tab').click();
        });
        
        templatesGrid.appendChild(templateCard);
    });
}

/**
 * 加载分子模板
 * @param {string} smiles - SMILES字符串
 */
function loadMoleculeTemplate(smiles) {
    if (window.jsmeApplet) {
        window.jsmeApplet.readGenericMolecularInput(smiles);
    }
}

/**
 * 监听深色模式变化
 */
function observeDarkModeChanges() {
    const observer = new MutationObserver(mutations => {
        mutations.forEach(mutation => {
            if (mutation.attributeName === 'data-theme') {
                updateEditorTheme();
            }
        });
    });

    observer.observe(document.documentElement, {
        attributes: true,
        attributeFilter: ['data-theme']
    });
}

/**
 * 更新编辑器主题
 */
function updateEditorTheme() {
    if (!window.jsmeApplet) return;

    const isDarkMode = document.documentElement.getAttribute('data-theme') === 'dark';
    const options = {
        ...JSME_OPTIONS,
        "borderColor": getComputedStyle(document.documentElement).getPropertyValue('--border-color'),
        "backgroundColor": getComputedStyle(document.documentElement).getPropertyValue('--bg-color')
    };

    if (isDarkMode) {
        options.atomColorPalette = "whiteOnBlack";
        options.bondColorPalette = "whiteOnBlack";
    }

    window.jsmeApplet.options(options);
    window.jsmeApplet.repaint();
}

/**
 * 清空编辑器
 */
window.clearEditor = function() {
    if (window.jsmeApplet) {
        window.jsmeApplet.reset();
    }
}

/**
 * 撤销上一步操作
 */
window.undoLastAction = function() {
    if (window.jsmeApplet) {
        window.jsmeApplet.undo();
    }
}

/**
 * 重做上一步操作
 */
window.redoLastAction = function() {
    if (window.jsmeApplet) {
        window.jsmeApplet.redo();
    }
}

/**
 * 从编辑器获取SMILES并计算
 */
window.calculateFromEditor = async function() {
    if (!window.jsmeApplet) {
        window.showToast('分子编辑器未加载完成，请稍候', 'warning');
        return;
    }

    const smiles = window.jsmeApplet.smiles();
    if (!smiles) {
        window.showToast('请先绘制分子结构', 'warning');
        return;
    }
    
    // 显示加载动画
    document.querySelector('.loading').style.display = 'block';
    document.getElementById('results').innerHTML = '';
    
    try {
        // 创建FormData
        const formData = new FormData();
        formData.append('smiles', smiles);
        
        // 发送请求
        const response = await fetch('/calculate', {
            method: 'POST',
            body: formData
        });
        
        if (!response.ok) {
            throw new Error('计算请求失败');
        }

        const data = await response.json();
        handleCalculationResponse(data);
    } catch (error) {
        window.showToast(error.message, 'error');
    } finally {
        document.querySelector('.loading').style.display = 'none';
    }
}

// 初始化SMILES示例点击事件
document.addEventListener('DOMContentLoaded', () => {
    document.querySelectorAll('.smiles-example').forEach(button => {
        button.addEventListener('click', () => {
            const smiles = button.getAttribute('data-smiles');
            document.getElementById('smiles').value = smiles;
        });
    });
});