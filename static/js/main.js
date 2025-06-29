// 主要JavaScript功能
class PiHOMACalculator {
    constructor() {
        this.history = JSON.parse(localStorage.getItem('piHomaHistory') || '[]');
        this.currentTheme = localStorage.getItem('theme') || 'light';
        this.init();
    }

    init() {
        this.initTheme();
        this.initEventListeners();
        this.initTooltips();
        this.loadHistory();
        this.initTemplates();
    }

    initTheme() {
        if (this.currentTheme === 'dark') {
            document.documentElement.setAttribute('data-theme', 'dark');
            const icon = document.querySelector('#darkModeToggle i');
            if (icon) {
                icon.className = 'fas fa-sun';
            }
        }
    }

    initEventListeners() {
        // 暗色模式切换
        document.getElementById('darkModeToggle')?.addEventListener('click', () => {
            this.toggleTheme();
        });

        // SMILES示例点击
        document.querySelectorAll('.smiles-example').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const smiles = e.target.dataset.smiles;
                document.getElementById('smiles').value = smiles;
            });
        });

        // 表单提交
        document.getElementById('calculate-form')?.addEventListener('submit', (e) => {
            e.preventDefault();
            this.calculateFromSmiles();
        });

        // 键盘快捷键
        document.addEventListener('keydown', (e) => {
            if (e.ctrlKey && e.key === 'Enter') {
                this.calculateFromActiveTab();
            }
        });
    }

    initTooltips() {
        // 初始化Bootstrap工具提示
        const tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]'));
        tooltipTriggerList.map(tooltipTriggerEl => {
            return new bootstrap.Tooltip(tooltipTriggerEl);
        });
    }

    initTemplates() {
        const templates = [
            { name: '苯', smiles: 'c1ccccc1', description: '最简单的芳香环' },
            { name: '吡啶', smiles: 'n1ccccc1', description: '含氮芳香环' },
            { name: '噻吩', smiles: 's1cccc1', description: '含硫五元环' },
            { name: '呋喃', smiles: 'o1cccc1', description: '含氧五元环' },
            { name: '萘', smiles: 'c1ccc2ccccc2c1', description: '双苯环系统' },
            { name: '蒽', smiles: 'c1ccc2cc3ccccc3cc2c1', description: '三苯环系统' },
            { name: '吡咯', smiles: '[nH]1cccc1', description: '含氮五元环' },
            { name: '咪唑', smiles: 'n1ccnc1', description: '含两个氮的五元环' },
            { name: '吡嗪', smiles: 'n1ccncc1', description: '1,4-二氮苯' },
            { name: '嘧啶', smiles: 'n1cnccc1', description: '1,3-二氮苯' }
        ];

        const templatesGrid = document.querySelector('.templates-grid');
        if (templatesGrid) {
            templatesGrid.innerHTML = templates.map(template => `
                <div class="template-card" data-smiles="${template.smiles}">
                    <h6>${template.name}</h6>
                    <p class="small text-muted">${template.description}</p>
                    <code class="small">${template.smiles}</code>
                </div>
            `).join('');

            // 添加模板点击事件
            templatesGrid.addEventListener('click', (e) => {
                const card = e.target.closest('.template-card');
                if (card) {
                    const smiles = card.dataset.smiles;
                    this.loadTemplate(smiles);
                }
            });
        }
    }

    loadTemplate(smiles) {
        // 切换到SMILES输入标签
        const smilesTab = document.getElementById('smiles-tab');
        const smilesInput = document.getElementById('smiles');
        
        if (smilesTab && smilesInput) {
            smilesTab.click();
            smilesInput.value = smiles;
            
            // 添加动画效果
            smilesInput.classList.add('animate__animated', 'animate__pulse');
            setTimeout(() => {
                smilesInput.classList.remove('animate__animated', 'animate__pulse');
            }, 1000);
        }
    }

    toggleTheme() {
        const currentTheme = document.documentElement.getAttribute('data-theme');
        const newTheme = currentTheme === 'dark' ? 'light' : 'dark';
        
        document.documentElement.setAttribute('data-theme', newTheme);
        localStorage.setItem('theme', newTheme);
        
        const icon = document.querySelector('#darkModeToggle i');
        if (icon) {
            icon.className = newTheme === 'dark' ? 'fas fa-sun' : 'fas fa-moon';
        }
        
        this.currentTheme = newTheme;
    }

    calculateFromActiveTab() {
        const activeTab = document.querySelector('.nav-link.active');
        if (activeTab) {
            const tabId = activeTab.id;
            switch (tabId) {
                case 'editor-tab':
                    this.calculateFromEditor();
                    break;
                case 'smiles-tab':
                    this.calculateFromSmiles();
                    break;
                default:
                    break;
            }
        }
    }

    async calculateFromSmiles() {
        const smilesInput = document.getElementById('smiles');
        const smiles = smilesInput?.value.trim();
        
        if (!smiles) {
            this.showError('请输入SMILES字符串');
            return;
        }

        await this.performCalculation({ smiles });
    }

    async calculateFromEditor() {
        // 这里需要与JSME编辑器集成
        if (typeof jsmeApplet !== 'undefined') {
            const smiles = jsmeApplet.smiles();
            if (!smiles) {
                this.showError('请在编辑器中绘制分子结构');
                return;
            }
            await this.performCalculation({ smiles });
        } else {
            this.showError('分子编辑器未初始化');
        }
    }

    async performCalculation(data) {
        this.showLoading(true);
        
        try {
            const response = await fetch('/calculate', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(data)
            });

            const result = await response.json();
            
            if (result.success) {
                this.displayResults(result);
                this.addToHistory(data.smiles, result);
            } else {
                this.showError(result.error || '计算失败');
            }
        } catch (error) {
            console.error('计算错误:', error);
            this.showError('网络错误，请稍后重试');
        } finally {
            this.showLoading(false);
        }
    }

    displayResults(result) {
        const resultsContainer = document.getElementById('results');
        if (!resultsContainer) return;

        const resultHtml = `
            <div class="col-lg-10">
                <div class="result-card animate__animated animate__fadeInUp">
                    <div class="card-header">
                        <h5 class="mb-0">
                            <i class="fas fa-chart-line me-2"></i>
                            π-HOMA计算结果
                        </h5>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <h6>分子结构</h6>
                                <img src="data:image/png;base64,${result.molecule_image}" 
                                     class="molecule-image" alt="分子结构图">
                            </div>
                            <div class="col-md-6">
                                <h6>π-HOMA值</h6>
                                <div class="table-responsive">
                                    <table class="table table-hover">
                                        <thead>
                                            <tr>
                                                <th>环编号</th>
                                                <th>原子索引</th>
                                                <th>π-HOMA值</th>
                                                <th>芳香性评价</th>
                                            </tr>
                                        </thead>
                                        <tbody>
                                            ${result.rings.map(ring => `
                                                <tr>
                                                    <td>环 ${ring.ring_number}</td>
                                                    <td class="atom-indices">[${ring.atoms.join(', ')}]</td>
                                                    <td class="pi-homa-value">${ring.pi_homa}</td>
                                                    <td>${this.getAromaticityDescription(ring.pi_homa)}</td>
                                                </tr>
                                            `).join('')}
                                        </tbody>
                                    </table>
                                </div>
                                <div class="mt-3">
                                    <small class="text-muted">
                                        <strong>说明：</strong>π-HOMA值范围0-1，值越接近1表示芳香性越强
                                    </small>
                                </div>
                            </div>
                        </div>
                        <div class="row mt-3">
                            <div class="col-12">
                                <div class="d-flex gap-2">
                                    <button class="btn btn-outline-primary" onclick="calculator.exportResults()">
                                        <i class="fas fa-download me-1"></i>导出结果
                                    </button>
                                    <button class="btn btn-outline-secondary" onclick="calculator.shareResults()">
                                        <i class="fas fa-share me-1"></i>分享结果
                                    </button>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        `;

        resultsContainer.innerHTML = resultHtml;
        resultsContainer.scrollIntoView({ behavior: 'smooth' });
    }

    getAromaticityDescription(piHoma) {
        if (piHoma >= 0.8) return '<span class="badge bg-success">强芳香性</span>';
        if (piHoma >= 0.5) return '<span class="badge bg-warning">中等芳香性</span>';
        if (piHoma >= 0.2) return '<span class="badge bg-info">弱芳香性</span>';
        return '<span class="badge bg-secondary">非芳香性</span>';
    }

    addToHistory(smiles, result) {
        const historyItem = {
            id: Date.now(),
            timestamp: new Date().toISOString(),
            smiles: smiles,
            result: result,
            date: new Date().toLocaleDateString('zh-CN')
        };

        this.history.unshift(historyItem);
        
        // 限制历史记录数量
        if (this.history.length > 50) {
            this.history = this.history.slice(0, 50);
        }

        localStorage.setItem('piHomaHistory', JSON.stringify(this.history));
        this.updateHistoryDisplay();
    }

    loadHistory() {
        this.updateHistoryDisplay();
    }

    updateHistoryDisplay() {
        const historyTableBody = document.getElementById('historyTableBody');
        if (!historyTableBody) return;

        if (this.history.length === 0) {
            historyTableBody.innerHTML = `
                <tr>
                    <td colspan="4" class="text-center text-muted">暂无计算历史</td>
                </tr>
            `;
            return;
        }

        historyTableBody.innerHTML = this.history.map(item => `
            <tr>
                <td>${item.date}</td>
                <td><code class="small">${item.smiles}</code></td>
                <td>
                    ${item.result.rings.map(ring => 
                        `环${ring.ring_number}: ${ring.pi_homa}`
                    ).join('<br>')}
                </td>
                <td>
                    <button class="btn btn-sm btn-outline-primary" 
                            onclick="calculator.loadFromHistory('${item.smiles}')">
                        <i class="fas fa-redo"></i>
                    </button>
                    <button class="btn btn-sm btn-outline-danger" 
                            onclick="calculator.deleteHistoryItem(${item.id})">
                        <i class="fas fa-trash"></i>
                    </button>
                </td>
            </tr>
        `).join('');
    }

    loadFromHistory(smiles) {
        document.getElementById('smiles').value = smiles;
        document.getElementById('smiles-tab').click();
        
        // 关闭历史模态框
        const modal = bootstrap.Modal.getInstance(document.getElementById('historyModal'));
        if (modal) modal.hide();
    }

    deleteHistoryItem(id) {
        this.history = this.history.filter(item => item.id !== id);
        localStorage.setItem('piHomaHistory', JSON.stringify(this.history));
        this.updateHistoryDisplay();
    }

    clearHistory() {
        if (confirm('确定要清空所有历史记录吗？')) {
            this.history = [];
            localStorage.removeItem('piHomaHistory');
            this.updateHistoryDisplay();
        }
    }

    exportHistory() {
        if (this.history.length === 0) {
            this.showError('没有历史记录可导出');
            return;
        }

        const csvContent = this.generateCSV();
        this.downloadFile(csvContent, 'pi-homa-history.csv', 'text/csv');
    }

    exportResults() {
        const resultsContainer = document.getElementById('results');
        if (!resultsContainer.innerHTML.trim()) {
            this.showError('没有结果可导出');
            return;
        }

        // 这里可以实现PDF导出或其他格式
        this.showInfo('导出功能开发中...');
    }

    shareResults() {
        if (navigator.share) {
            navigator.share({
                title: 'π-HOMA计算结果',
                text: '查看我的分子芳香性计算结果',
                url: window.location.href
            });
        } else {
            // 复制链接到剪贴板
            navigator.clipboard.writeText(window.location.href).then(() => {
                this.showSuccess('链接已复制到剪贴板');
            });
        }
    }

    generateCSV() {
        const headers = ['日期', 'SMILES', '环编号', 'π-HOMA值'];
        const rows = [headers.join(',')];

        this.history.forEach(item => {
            item.result.rings.forEach(ring => {
                const row = [
                    item.date,
                    `"${item.smiles}"`,
                    ring.ring_number,
                    ring.pi_homa
                ];
                rows.push(row.join(','));
            });
        });

        return rows.join('\n');
    }

    downloadFile(content, filename, mimeType) {
        const blob = new Blob([content], { type: mimeType });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(url);
    }

    showLoading(show) {
        const loading = document.querySelector('.loading');
        if (loading) {
            loading.style.display = show ? 'block' : 'none';
        }
    }

    showError(message) {
        this.showToast(message, 'error');
    }

    showSuccess(message) {
        this.showToast(message, 'success');
    }

    showInfo(message) {
        this.showToast(message, 'info');
    }

    showToast(message, type = 'info') {
        // 创建toast容器（如果不存在）
        let toastContainer = document.getElementById('toast-container');
        if (!toastContainer) {
            toastContainer = document.createElement('div');
            toastContainer.id = 'toast-container';
            toastContainer.className = 'toast-container position-fixed top-0 end-0 p-3';
            toastContainer.style.zIndex = '1055';
            document.body.appendChild(toastContainer);
        }

        const toastId = 'toast-' + Date.now();
        const iconClass = {
            success: 'fas fa-check-circle text-success',
            error: 'fas fa-exclamation-circle text-danger',
            info: 'fas fa-info-circle text-info',
            warning: 'fas fa-exclamation-triangle text-warning'
        }[type] || 'fas fa-info-circle text-info';

        const toastHtml = `
            <div id="${toastId}" class="toast" role="alert">
                <div class="toast-header">
                    <i class="${iconClass} me-2"></i>
                    <strong class="me-auto">π-HOMA计算器</strong>
                    <button type="button" class="btn-close" data-bs-dismiss="toast"></button>
                </div>
                <div class="toast-body">
                    ${message}
                </div>
            </div>
        `;

        toastContainer.insertAdjacentHTML('beforeend', toastHtml);
        
        const toastElement = document.getElementById(toastId);
        const toast = new bootstrap.Toast(toastElement, {
            autohide: true,
            delay: type === 'error' ? 5000 : 3000
        });
        
        toast.show();

        // 清理已隐藏的toast
        toastElement.addEventListener('hidden.bs.toast', () => {
            toastElement.remove();
        });
    }
}

// 全局函数（保持向后兼容）
function calculateFromSmiles() {
    calculator.calculateFromSmiles();
}

function calculateFromEditor() {
    calculator.calculateFromEditor();
}

function clearHistory() {
    calculator.clearHistory();
}

function exportHistory() {
    calculator.exportHistory();
}

// 初始化应用
let calculator;
document.addEventListener('DOMContentLoaded', () => {
    calculator = new PiHOMACalculator();
});