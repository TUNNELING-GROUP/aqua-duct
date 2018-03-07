"nnoremap <F9> :silent !ctags-exuberant $(find src -iname '*.py')<CR> :redraw!<CR>
nnoremap <F9> :silent !ectags $(find src -iname '*.py')<CR> :redraw!<CR>
set t_Co=256
colorscheme dracula
"autocmd VimEnter * NERDTreeTabsOpen
nnoremap <F12> :silent !bash -c './make_pdf.sh &>/dev/null' &<CR> :redraw!<CR>

