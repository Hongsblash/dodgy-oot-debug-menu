glabel func_80B9DF1C
/* 00D0C 80B9DF1C 908E0191 */  lbu     $t6, 0x0191($a0)           ## 00000191
/* 00D10 80B9DF20 31CF0002 */  andi    $t7, $t6, 0x0002           ## $t7 = 00000000
/* 00D14 80B9DF24 51E00017 */  beql    $t7, $zero, .L80B9DF84     
/* 00D18 80B9DF28 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 00D1C 80B9DF2C 9098017F */  lbu     $t8, 0x017F($a0)           ## 0000017F
/* 00D20 80B9DF30 33190002 */  andi    $t9, $t8, 0x0002           ## $t9 = 00000000
/* 00D24 80B9DF34 57200013 */  bnel    $t9, $zero, .L80B9DF84     
/* 00D28 80B9DF38 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 00D2C 80B9DF3C 8C820188 */  lw      $v0, 0x0188($a0)           ## 00000188
/* 00D30 80B9DF40 50400010 */  beql    $v0, $zero, .L80B9DF84     
/* 00D34 80B9DF44 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 00D38 80B9DF48 84480032 */  lh      $t0, 0x0032($v0)           ## 00000032
/* 00D3C 80B9DF4C 848900B6 */  lh      $t1, 0x00B6($a0)           ## 000000B6
/* 00D40 80B9DF50 01091823 */  subu    $v1, $t0, $t1              
/* 00D44 80B9DF54 00031C00 */  sll     $v1, $v1, 16               
/* 00D48 80B9DF58 00031C03 */  sra     $v1, $v1, 16               
/* 00D4C 80B9DF5C 04600003 */  bltz    $v1, .L80B9DF6C            
/* 00D50 80B9DF60 00031023 */  subu    $v0, $zero, $v1            
/* 00D54 80B9DF64 10000001 */  beq     $zero, $zero, .L80B9DF6C   
/* 00D58 80B9DF68 00601025 */  or      $v0, $v1, $zero            ## $v0 = 00000000
.L80B9DF6C:
/* 00D5C 80B9DF6C 28415001 */  slti    $at, $v0, 0x5001           
/* 00D60 80B9DF70 54200004 */  bnel    $at, $zero, .L80B9DF84     
/* 00D64 80B9DF74 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 00D68 80B9DF78 03E00008 */  jr      $ra                        
/* 00D6C 80B9DF7C 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L80B9DF80:
/* 00D70 80B9DF80 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L80B9DF84:
/* 00D74 80B9DF84 03E00008 */  jr      $ra                        
/* 00D78 80B9DF88 00000000 */  nop


