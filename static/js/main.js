/**
 * 主要的应用逻辑
 */

// 当文档加载完成时初始化
document.addEventListener('DOMContentLoaded', function() {
    initializeFormHandler();
    initializeDarkMode();
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