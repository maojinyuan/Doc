-- Options are automatically loaded before lazy.nvim startup
-- Default options that are always set: https://github.com/LazyVim/LazyVim/blob/main/lua/lazyvim/config/options.lua
-- Add any additional options here

local opt = vim.opt
opt.relativenumber = false
opt.mouse="c"
opt.undofile = false
opt.virtualedit = ""
-- vim.api.nvim_set_keymap('i', '<A-e>', '<M-w><CR>', { noremap = true, silent = true })
vim.api.nvim_set_keymap('n', '<RightMouse>', '<NOP>', { noremap = true, silent = true })
-- 自定义函数，用于在新行插入文本时不继承上一行的引号
function OpenNewLine()
    local col = vim.fn.col('.')  -- 获取当前光标列位置
    local line = vim.fn.line('.')  -- 获取当前光标行位置

    -- 在新行插入一个换行符
    vim.fn.append(line, '')

    -- 移动到新行的起始位置
    vim.fn.cursor(line + 1, col)

    -- 如果当前行有引号，则删除引号
    local current_line = vim.fn.getline(line)
    if string.find(current_line, '[\'"]', 1) then
        vim.fn.setline(line, string.gsub(current_line, '[\'"]', ''))
    end

    -- 切换到插入模式
    vim.fn.feedkeys('i')
end
vim.api.nvim_set_keymap('n', 'o', ':lua OpenNewLine()<CR>', { noremap = true, silent = true })

