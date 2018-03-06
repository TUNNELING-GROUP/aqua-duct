nnoremap <F9> :silent !ctags-exuberant $(find src -iname '*.py')<CR> :redraw!<CR>
set t_Co=256
colorscheme dracula
nnoremap <F12> :silent !./docs/make_pdf.sh &> /dev/null&<CR> :redraw!<CR>

