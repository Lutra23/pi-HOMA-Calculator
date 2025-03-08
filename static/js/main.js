/**
 * 主要的应用逻辑
 */

// 全局变量
let darkMode = false;
let calculationHistory = [];

// DOM加载完成后初始化
document.addEventListener('DOMContentLoaded', () => {
    console.log('DOM加载完成，开始初始化...');
    
    // 初始化工具提示
    const tooltips = document.querySelectorAll('[data-bs-toggle="tooltip"]');
    tooltips.forEach(tooltip => new bootstrap.Tooltip(tooltip));

    // 初始化暗色模式
    initDarkMode();

    // 初始化历史记录
    loadHistory();

    // 绑定事件监听器
    bindEventListeners();
});

/**
 * 初始化深色模式
 */
function initializeDarkMode() {
    const darkModeToggle = document.getElementById('darkModeToggle');
    const prefersDarkScheme = window.matchMedia('(prefers-color-scheme: dark)');
    
    // 检查本地存储中的主题设置
    const currentTheme = localStorage.getItem('theme');
    if (currentTheme) {
        document.documentElement.setAttribute('data-theme', currentTheme);
        updateDarkModeIcon(currentTheme === 'dark');
    } else if (prefersDarkScheme.matches) {
        document.documentElement.setAttribute('data-theme', 'dark');
        updateDarkModeIcon(true);
    }
    
    // 监听系统主题变化
    prefersDarkScheme.addEventListener('change', (e) => {
        if (!localStorage.getItem('theme')) {
            const newTheme = e.matches ? 'dark' : 'light';
            document.documentElement.setAttribute('data-theme', newTheme);
            updateDarkModeIcon(e.matches);
        }
    });
    
    // 切换按钮点击事件
    darkModeToggle.addEventListener('click', () => {
        const currentTheme = document.documentElement.getAttribute('data-theme');
        const newTheme = currentTheme === 'dark' ? 'light' : 'dark';
        
        document.documentElement.setAttribute('data-theme', newTheme);
        localStorage.setItem('theme', newTheme);
        updateDarkModeIcon(newTheme === 'dark');
    });
}

/**
 * 更新深色模式图标
 * @param {boolean} isDark
 */
function updateDarkModeIcon(isDark) {
    const icon = document.querySelector('#darkModeToggle i');
    icon.className = isDark ? 'fas fa-sun' : 'fas fa-moon';
}

/**
 * 初始化表单处理器
 */
function initializeFormHandler() {
    const form = document.getElementById('calculate-form');
    if (!form) return;

    form.addEventListener('submit', async function(e) {
        e.preventDefault();
        
        const smilesInput = document.getElementById('smiles');
        if (!smilesInput.value.trim()) {
            window.showToast('请输入SMILES字符串', 'warning');
            return;
        }
        
        // 显示加载动画
        document.querySelector('.loading').style.display = 'block';
        document.getElementById('results').innerHTML = '';
        
        try {
            // 发送请求
            const response = await fetch('/calculate', {
                method: 'POST',
                body: new FormData(this)
            });
            
            if (!response.ok) {
                throw new Error('计算请求失败');
            }

            const data = await response.json();
            handleCalculationResponse(data);
        } catch (error) {
            window.showToast(error.message, 'error');
            showError(error.message);
        } finally {
            document.querySelector('.loading').style.display = 'none';
        }
    });
}

/**
 * 处理计算响应
 * @param {Object} data - 服务器返回的响应数据
 */
window.handleCalculationResponse = function(data) {
    if (data.success) {
        const resultsHtml = `
            <div class="col-md-8 fade-in">
                <div class="card result-card">
                    <div class="card-body">
                        <h5 class="card-title">
                            <i class="fas fa-chart-bar me-2"></i>计算结果
                            <button class="btn btn-sm btn-outline-primary float-end" onclick="exportResults(this)">
                                <i class="fas fa-download me-1"></i>导出
                            </button>
                        </h5>
                        <div class="text-center">
                            <img src="data:image/png;base64,${data.molecule_image}" 
                                 class="img-fluid molecule-image" alt="分子结构">
                            <div class="legend">
                                <i class="fas fa-info-circle me-1"></i>
                                <small>不同颜色表示不同的环系统，数字为原子编号</small>
                            </div>
                        </div>
                        <div class="ring-info mt-4">
                            <h6 class="mb-3">
                                <i class="fas fa-ring me-2"></i>环系统分析
                            </h6>
                            <div class="table-responsive">
                                <table class="table table-hover">
                                    <thead>
                                        <tr>
                                            <th>环编号</th>
                                            <th>原子索引</th>
                                            <th>π-HOMA值</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        ${data.rings.map(ring => `
                                            <tr>
                                                <td>${ring.ring_number}</td>
                                                <td class="atom-indices">${ring.atoms.join(', ')}</td>
                                                <td class="pi-homa-value">${ring.pi_homa}</td>
                                            </tr>
                                        `).join('')}
                                    </tbody>
                                </table>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        `;
        
        document.getElementById('results').innerHTML = resultsHtml;
        window.showToast('计算完成！', 'success');
    } else {
        showError(data.error);
    }
}

/**
 * 显示错误信息
 * @param {string} error - 错误信息
 */
function showError(error) {
    document.getElementById('results').innerHTML = `
        <div class="col-md-8 fade-in">
            <div class="alert alert-danger d-flex align-items-center" role="alert">
                <i class="fas fa-exclamation-circle me-2"></i>
                <div>计算过程中出错：${error}</div>
            </div>
        </div>
    `;
}

/**
 * 导出计算结果
 * @param {HTMLElement} button - 导出按钮元素
 */
window.exportResults = function(button) {
    const resultCard = button.closest('.result-card');
    if (!resultCard) return;

    // 临时禁用按钮
    button.disabled = true;
    const originalText = button.innerHTML;
    button.innerHTML = '<i class="fas fa-spinner fa-spin me-1"></i>导出中...';

    import('html2canvas').then(html2canvas => {
        html2canvas.default(resultCard).then(canvas => {
            // 转换为图片
            const imgData = canvas.toDataURL('image/png');
            
            // 创建下载链接
            const link = document.createElement('a');
            link.download = `pi-homa-results-${new Date().getTime()}.png`;
            link.href = imgData;
            link.click();

            // 恢复按钮状态
            button.disabled = false;
            button.innerHTML = originalText;
            
            window.showToast('导出成功！', 'success');
        });
    }).catch(error => {
        button.disabled = false;
        button.innerHTML = originalText;
        window.showToast('导出失败，请重试', 'error');
        console.error('导出错误：', error);
    });
}

// 初始化暗色模式
function initDarkMode() {
    // 从localStorage读取暗色模式设置
    darkMode = localStorage.getItem('darkMode') === 'true';
    updateDarkMode();
}

// 更新暗色模式
function updateDarkMode() {
    document.documentElement.setAttribute('data-theme', darkMode ? 'dark' : 'light');
    const icon = document.querySelector('#darkModeToggle i');
    icon.className = darkMode ? 'fas fa-sun' : 'fas fa-moon';
    localStorage.setItem('darkMode', darkMode);
}

// 绑定事件监听器
function bindEventListeners() {
    console.log('开始绑定事件监听器...');
    
    // 暗色模式切换
    const darkModeToggle = document.getElementById('darkModeToggle');
    if (darkModeToggle) {
        darkModeToggle.addEventListener('click', () => {
            darkMode = !darkMode;
            updateDarkMode();
        });
    }

    // SMILES示例点击
    document.querySelectorAll('.smiles-example').forEach(button => {
        button.addEventListener('click', () => {
            const smiles = button.dataset.smiles;
            const smilesInput = document.getElementById('smiles');
            if (smilesInput) {
                smilesInput.value = smiles;
                console.log('设置SMILES示例:', smiles);
            }
        });
    });

    // 计算按钮点击事件
    const calculateBtn = document.getElementById('calculate-btn');
    if (calculateBtn) {
        console.log('找到计算按钮，绑定点击事件');
        calculateBtn.onclick = async (e) => {
            e.preventDefault();
            console.log('计算按钮被点击');
            await calculateFromSmiles();
        };
    } else {
        console.error('未找到计算按钮');
    }
}

// 从SMILES计算
window.calculateFromSmiles = async function() {
    console.log('开始SMILES计算');
    const smilesInput = document.getElementById('smiles');
    if (!smilesInput) {
        console.error('未找到SMILES输入框');
        return;
    }
    
    const smiles = smilesInput.value.trim();
    console.log('SMILES输入值:', smiles);
    
    if (!smiles) {
        showToast('错误', '请输入SMILES字符串');
        return;
    }

    showLoading();
    try {
        console.log('发送计算请求...');
        const formData = new FormData();
        formData.append('smiles', smiles);
        
        const response = await fetch('/calculate', {
            method: 'POST',
            body: formData
        });

        console.log('收到服务器响应:', response.status);
        if (!response.ok) {
            throw new Error('计算请求失败');
        }

        const result = await response.json();
        console.log('解析响应数据:', result);
        
        if (!result.success) {
            throw new Error(result.error || '计算失败');
        }
        
        if (!result.rings || !Array.isArray(result.rings)) {
            throw new Error('返回数据格式错误: 缺少环系统数据');
        }
        
        displayResult(result);
        addToHistory(result);
        showToast('成功', '计算完成！');
    } catch (error) {
        console.error('计算错误:', error);
        showToast('错误', error.message);
        const resultsDiv = document.getElementById('results');
        if (resultsDiv) {
            resultsDiv.innerHTML = `
                <div class="col-lg-10">
                    <div class="alert alert-danger">
                        <h5><i class="fas fa-exclamation-triangle me-2"></i>计算出错</h5>
                        <p>${error.message}</p>
                    </div>
                </div>
            `;
        }
    } finally {
        hideLoading();
    }
}

// 显示结果
function displayResult(result) {
    console.log('开始显示结果，数据:', result);
    
    if (!result.rings || !Array.isArray(result.rings) || result.rings.length === 0) {
        throw new Error('没有找到芳香环系统');
    }
    
    const resultsDiv = document.getElementById('results');
    if (!resultsDiv) {
        console.error('未找到结果显示区域');
        return;
    }
    
    // 处理环系统数据
    const processedRings = result.rings.map(ring => {
        // 确保所有必要的字段都存在
        const ringNumber = ring.ring_number || '未知';
        const atoms = Array.isArray(ring.atoms) ? ring.atoms : [];
        // 处理pi_homa值，确保是数字并保留4位小数
        let piHoma = '未知';
        if (ring.pi_homa !== undefined && ring.pi_homa !== null) {
            const piHomaNum = parseFloat(ring.pi_homa);
            if (!isNaN(piHomaNum)) {
                piHoma = piHomaNum.toFixed(4);
            }
        }
        return { ringNumber, atoms, piHoma };
    });
    
    console.log('处理后的环系统数据:', processedRings);
    
    const resultHtml = `
        <div class="col-lg-10">
            <div class="result-card fade-in-up">
                <div class="card-header">
                    <h5 class="mb-0">计算结果</h5>
                </div>
                <div class="card-body">
                    <div class="row">
                        <div class="col-md-6">
                            <h6>分子结构</h6>
                            <img src="data:image/png;base64,${result.molecule_image}" 
                                 alt="分子结构" class="molecule-image">
                        </div>
                        <div class="col-md-6">
                            <h6>环系统分析</h6>
                            <div class="table-responsive">
                                <table class="table table-hover">
                                    <thead>
                                        <tr>
                                            <th>环编号</th>
                                            <th>原子索引</th>
                                            <th>π-HOMA值</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        ${processedRings.map(ring => `
                                            <tr>
                                                <td>${ring.ringNumber}</td>
                                                <td class="atom-indices">${ring.atoms.join(', ')}</td>
                                                <td class="pi-homa-value">${ring.piHoma}</td>
                                            </tr>
                                        `).join('')}
                                    </tbody>
                                </table>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    `;
    
    console.log('更新结果显示');
    resultsDiv.innerHTML = resultHtml;
}

// 添加到历史记录
function addToHistory(result) {
    const timestamp = new Date().toLocaleString();
    const historyItem = {
        timestamp,
        smiles: result.smiles,
        rings: result.rings.map(ring => ({
            ring_number: ring.ring_number,
            atoms: ring.atoms,
            pi_homa: ring.pi_homa
        })),
        molecule_image: result.molecule_image
    };
    
    calculationHistory.unshift(historyItem);
    if (calculationHistory.length > 50) {
        calculationHistory.pop();
    }
    
    saveHistory();
    updateHistoryTable();
}

// 保存历史记录
function saveHistory() {
    localStorage.setItem('calculationHistory', JSON.stringify(calculationHistory));
}

// 加载历史记录
function loadHistory() {
    const saved = localStorage.getItem('calculationHistory');
    if (saved) {
        calculationHistory = JSON.parse(saved);
        updateHistoryTable();
    }
}

// 更新历史记录表格
function updateHistoryTable() {
    const tbody = document.getElementById('historyTableBody');
    if (!tbody) return;

    tbody.innerHTML = calculationHistory.map((item, index) => {
        // 获取第一个环的π-HOMA值作为显示（如果有多个环，只显示第一个）
        const firstRing = item.rings && item.rings[0];
        const piHoma = firstRing && firstRing.pi_homa !== undefined ? 
            parseFloat(firstRing.pi_homa).toFixed(4) : '未知';
        
        return `
            <tr>
                <td>${item.timestamp}</td>
                <td>
                    <img src="data:image/png;base64,${item.molecule_image}" 
                         alt="分子结构" style="height: 50px;">
                </td>
                <td class="pi-homa-value">${piHoma}</td>
                <td>
                    <button class="btn btn-sm btn-outline-primary" 
                            onclick="loadHistoryItem(${index})">
                        <i class="fas fa-redo"></i>
                    </button>
                </td>
            </tr>
        `;
    }).join('');
}

// 加载历史记录项
function loadHistoryItem(index) {
    const item = calculationHistory[index];
    document.getElementById('smiles').value = item.smiles;
    document.getElementById('smiles-tab').click();
}

// 清空历史记录
function clearHistory() {
    if (confirm('确定要清空所有历史记录吗？')) {
        calculationHistory = [];
        saveHistory();
        updateHistoryTable();
    }
}

// 导出历史记录
function exportHistory() {
    const csv = [
        ['时间', 'SMILES', '环系统数', 'π-HOMA值(所有环)'].join(','),
        ...calculationHistory.map(item => [
            item.timestamp,
            item.smiles,
            item.rings ? item.rings.length : 0,
            item.rings ? item.rings.map(r => r.pi_homa).join(';') : ''
        ].join(','))
    ].join('\n');

    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = 'pi-homa-history.csv';
    link.click();
}

// 导出分子结构
function exportMolecule() {
    const smiles = jsmeApplet.smiles();
    if (!smiles) {
        showToast('错误', '请先绘制分子结构');
        return;
    }

    const blob = new Blob([smiles], { type: 'text/plain;charset=utf-8;' });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = 'molecule.smi';
    link.click();
}

// 显示加载动画
function showLoading() {
    document.querySelector('.loading').style.display = 'block';
}

// 隐藏加载动画
function hideLoading() {
    document.querySelector('.loading').style.display = 'none';
}

// 显示提示消息
function showToast(title, message) {
    console.log('显示提示:', title, message);
    const toastHtml = `
        <div class="toast-container position-fixed bottom-0 end-0 p-3">
            <div class="toast" role="alert">
                <div class="toast-header">
                    <strong class="me-auto">${title}</strong>
                    <button type="button" class="btn-close" data-bs-dismiss="toast"></button>
                </div>
                <div class="toast-body">${message}</div>
            </div>
        </div>
    `;

    document.body.insertAdjacentHTML('beforeend', toastHtml);
    const toastEl = document.querySelector('.toast');
    const toast = new bootstrap.Toast(toastEl);
    toast.show();

    toastEl.addEventListener('hidden.bs.toast', () => {
        toastEl.parentElement.remove();
    });
}