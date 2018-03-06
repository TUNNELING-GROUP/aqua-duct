nnoremap <F9> :silent !ctags-exuberant $(find src -iname '*.py')<CR> :redraw!<CR>
set t_Co=256
colorscheme dracula
autocmd VimEnter * NERDTreeTabsOpen

