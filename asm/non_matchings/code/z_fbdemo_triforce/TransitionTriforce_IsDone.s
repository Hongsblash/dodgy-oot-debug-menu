.late_rodata
glabel D_8014309C
    .float 0.03

.text
glabel TransitionTriforce_IsDone
/* B29AE4 800B2944 8C82000C */  lw    $v0, 0xc($a0)
/* B29AE8 800B2948 24010001 */  li    $at, 1
/* B29AEC 800B294C 10410002 */  beq   $v0, $at, .L800B2958
/* B29AF0 800B2950 24010002 */   li    $at, 2
/* B29AF4 800B2954 1441000B */  bne   $v0, $at, .L800B2984
.L800B2958:
/* B29AF8 800B2958 3C018014 */   lui   $at, %hi(D_8014309C)
/* B29AFC 800B295C C424309C */  lwc1  $f4, %lo(D_8014309C)($at)
/* B29B00 800B2960 C4860004 */  lwc1  $f6, 4($a0)
/* B29B04 800B2964 00001025 */  move  $v0, $zero
/* B29B08 800B2968 4604303E */  c.le.s $f6, $f4
/* B29B0C 800B296C 00000000 */  nop
/* B29B10 800B2970 45000002 */  bc1f  .L800B297C
/* B29B14 800B2974 00000000 */   nop
/* B29B18 800B2978 24020001 */  li    $v0, 1
.L800B297C:
/* B29B1C 800B297C 03E00008 */  jr    $ra
/* B29B20 800B2980 00000000 */   nop

.L800B2984:
/* B29B24 800B2984 24010003 */  li    $at, 3
/* B29B28 800B2988 10410002 */  beq   $v0, $at, .L800B2994
/* B29B2C 800B298C 24010004 */   li    $at, 4
/* B29B30 800B2990 1441000B */  bne   $v0, $at, .L800B29C0
.L800B2994:
/* B29B34 800B2994 3C013F80 */   li    $at, 0x3F800000 # 0.000000
/* B29B38 800B2998 44815000 */  mtc1  $at, $f10
/* B29B3C 800B299C C4880004 */  lwc1  $f8, 4($a0)
/* B29B40 800B29A0 00001025 */  move  $v0, $zero
/* B29B44 800B29A4 4608503E */  c.le.s $f10, $f8
/* B29B48 800B29A8 00000000 */  nop
/* B29B4C 800B29AC 45000002 */  bc1f  .L800B29B8
/* B29B50 800B29B0 00000000 */   nop
/* B29B54 800B29B4 24020001 */  li    $v0, 1
.L800B29B8:
/* B29B58 800B29B8 03E00008 */  jr    $ra
/* B29B5C 800B29BC 00000000 */   nop

.L800B29C0:
/* B29B60 800B29C0 00001025 */  move  $v0, $zero
/* B29B64 800B29C4 03E00008 */  jr    $ra
/* B29B68 800B29C8 00000000 */   nop
