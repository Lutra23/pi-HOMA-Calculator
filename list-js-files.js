const fs = require('fs');
const path = require('path');

function findJsFiles(dir, maxResults = 5) {
  const results = [];
  
  function searchDirectory(currentDir) {
    if (results.length >= maxResults) return;
    
    try {
      const items = fs.readdirSync(currentDir);
      
      for (const item of items) {
        if (results.length >= maxResults) break;
        
        const fullPath = path.join(currentDir, item);
        const stat = fs.statSync(fullPath);
        
        if (stat.isDirectory()) {
          searchDirectory(fullPath);
        } else if (stat.isFile() && item.endsWith('.js')) {
          results.push(fullPath);
        }
      }
    } catch (error) {
      // Skip directories that can't be read
      console.error(`Error reading directory ${currentDir}:`, error.message);
    }
  }
  
  searchDirectory(dir);
  return results;
}

// Find JavaScript files in static/jsme directory
const jsmeDir = 'static/jsme';
const jsFiles = findJsFiles(jsmeDir, 5);

console.log('JavaScript files found in static/jsme:');
jsFiles.forEach(file => console.log(file));

if (jsFiles.length === 0) {
  console.log('No JavaScript files found in the static/jsme directory.');
}