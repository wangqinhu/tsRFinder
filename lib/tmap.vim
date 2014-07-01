" Vim syntax file
" Language:     tsRFinder tmap
" Maintainer:   Qinhu Wang <wangqinhu@nwafu.edu.cn>
" Last Change:  Jue 30, 2014

" For version 5.x: Clear all syntax items
" For version 6.x: Quit when a syntax file was already loaded
if version < 600
    syntax clear
elseif exists("b:current_syntax")
    finish
endif

" tmap syn 
syn match baseA       "[Aa]"
syn match baseT       "[Tt]"
syn match baseC       "[Cc]"
syn match baseG       "[Gg]"
syn match gap         "-"
syn match clover      "[<.>]"
syn match bdi         "\d"
syn match head        "<tmap"
syn match tail        "/>"
syn region idnum start=/\s\+/ end=/\S\+$/

" tmap highlight
highlight baseA         ctermfg=Green          guifg=Green
highlight baseT         ctermfg=Red            guifg=Red
highlight baseC         ctermfg=Blue           guifg=Blue
highlight baseG         ctermfg=Yellow         guifg=Yellow
highlight gap           ctermfg=Grey           guifg=Grey
highlight clover        ctermfg=LightGreen     guifg=DarkYellow
highlight bdi           ctermfg=LightRed       guifg=LightRed
highlight head          ctermfg=DarkYellow     guifg=DarkYellow
highlight tail          ctermfg=DarkYellow     guifg=DarkYellow
highlight idnum         ctermfg=LightBlue      guifg=DarkGreen

let b:current_syntax = "tmap"
