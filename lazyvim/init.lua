-- bootstrap lazy.nvim, LazyVim and your plugins
require("config.lazy")
----------

-- 禁用自动补全
local cmp = require('cmp')

-- Configure nvim-cmp
cmp.setup {
  -- other configuration options

  -- Disable auto-completion
  completion = {
    autocomplete = false
  }
}

