glabel func_80AAD518
/* 00C18 80AAD518 2401000B */  addiu   $at, $zero, 0x000B         ## $at = 0000000B
/* 00C1C 80AAD51C AFA40000 */  sw      $a0, 0x0000($sp)           
/* 00C20 80AAD520 AFA60008 */  sw      $a2, 0x0008($sp)           
/* 00C24 80AAD524 14A1000B */  bne     $a1, $at, .L80AAD554       
/* 00C28 80AAD528 AFA7000C */  sw      $a3, 0x000C($sp)           
/* 00C2C 80AAD52C 8FA20014 */  lw      $v0, 0x0014($sp)           
/* 00C30 80AAD530 8FA30010 */  lw      $v1, 0x0010($sp)           
/* 00C34 80AAD534 844F027A */  lh      $t7, 0x027A($v0)           ## 0000027A
/* 00C38 80AAD538 846E0002 */  lh      $t6, 0x0002($v1)           ## 00000002
/* 00C3C 80AAD53C 84790004 */  lh      $t9, 0x0004($v1)           ## 00000004
/* 00C40 80AAD540 01CFC023 */  subu    $t8, $t6, $t7              
/* 00C44 80AAD544 A4780002 */  sh      $t8, 0x0002($v1)           ## 00000002
/* 00C48 80AAD548 84480278 */  lh      $t0, 0x0278($v0)           ## 00000278
/* 00C4C 80AAD54C 03284821 */  addu    $t1, $t9, $t0              
/* 00C50 80AAD550 A4690004 */  sh      $t1, 0x0004($v1)           ## 00000004
.L80AAD554:
/* 00C54 80AAD554 03E00008 */  jr      $ra                        
/* 00C58 80AAD558 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
