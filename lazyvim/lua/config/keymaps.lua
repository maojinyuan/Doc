-- Keymaps are automatically loaded on the VeryLazy event
-- Default keymaps that are always set: https://github.com/LazyVim/LazyVim/blob/main/lua/lazyvim/config/keymaps.lua
-- Add any additional keymaps here
--
local keymap = vim.keymap
keymap.set('n', 'wq', ':wq <CR>')
-- keymap.set('n', '<leader>s', ':w <CR>')
keymap.set('n', '<leader>j', '<c-w>j')
keymap.set('n', '<leader>k', '<c-w>k')
keymap.set('n', '<leader>h', '<c-w>h')
keymap.set('n', '<leader>l', '<c-w>l')

