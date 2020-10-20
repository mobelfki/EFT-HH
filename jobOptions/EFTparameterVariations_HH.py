def varyEFTParam(par, val):

    print "Vary parameter nb. ", par

    non_zero_params = dict()

    if par == 0:
        print "Compute SM-limit (all EFT param set to 0)."

    elif par == 2:
        non_zero_params={'cpDC':val}
    elif par == 3:
        non_zero_params={'cpWB':val}
    elif par == 4:
        non_zero_params={'cdp':val}
    elif par == 5:
        non_zero_params={'cp':val}
    elif par == 6:
        non_zero_params={'cWWW':val}  # scaling to Higgs cross-section 
    elif par == 7:
        non_zero_params={'cG':val}
    elif par == 8:
        non_zero_params={'cpG':val}
    elif par == 9:
        non_zero_params={'cpW':val}
    elif par == 10:
        non_zero_params={'cpBB':val}
    elif par == 11:
        non_zero_params={'cHW':val}
    elif par == 12:
        non_zero_params={'cHWtil':val}
    elif par == 13:
        non_zero_params={'cHB':val}
    elif par == 14:
        non_zero_params={'cHBtil':val}
    elif par == 15:
        non_zero_params={'cHWB':val}
    elif par == 16:
        non_zero_params={'cHWBtil':val}
    elif par == 17:
        non_zero_params={'ceHAbs':val}
    elif par == 18:
        non_zero_params={'cuHAbs':val}
    elif par == 19:
        non_zero_params={'cdHAbs':val}
    elif par == 20:
        non_zero_params={'cpqMi':val}
    elif par == 38:
        non_zero_params={'ceBAbs':val}
    elif par == 39:
        non_zero_params={'cuGAbs':val}
    elif par == 40:
        non_zero_params={'cuWAbs':val}
    elif par == 41:
        non_zero_params={'cuBAbs':val}
    elif par == 42:
        non_zero_params={'cdGAbs':val}
    elif par == 43:
        non_zero_params={'cdWAbs':val}
    elif par == 44:
        non_zero_params={'cdBAbs':val}
    elif par == 45:
        non_zero_params={'cHl1':val}
    elif par == 46:
        non_zero_params={'cHl3':val}
    elif par == 47:
        non_zero_params={'cHe':val}
    elif par == 48:
        non_zero_params={'cHq1':val}
    elif par == 49:
        non_zero_params={'cHq3':val}
    elif par == 50:
        non_zero_params={'cHu':val}
    elif par == 51:
        non_zero_params={'cHd':val}
    elif par == 52:
        non_zero_params={'cHudAbs':val}
    elif par == 53:
        non_zero_params={'cll':val}
    elif par == 54:
        non_zero_params={'cll1':val}
    elif par == 55:
        non_zero_params={'cqq1':val}
    elif par == 56:
        non_zero_params={'cqq11':val}
    elif par == 57:
        non_zero_params={'cqq3':val}
    elif par == 58:
        non_zero_params={'cqq31':val}
    elif par == 59:
        non_zero_params={'clq1':val}
    elif par == 60:
        non_zero_params={'clq3':val}
    elif par == 61:
        non_zero_params={'cee':val}
    elif par == 62:
        non_zero_params={'cuu':val}
    elif par == 63:
        non_zero_params={'cuu1':val}
    elif par == 64:
        non_zero_params={'cdd':val}
    elif par == 65:
        non_zero_params={'cdd1':val}
    elif par == 66:
        non_zero_params={'ceu':val}
    elif par == 67:
        non_zero_params={'ced':val}
    elif par == 68:
        non_zero_params={'cud1':val}
    elif par == 69:
        non_zero_params={'cud8':val}
    elif par == 70:
        non_zero_params={'cle':val}
    elif par == 71:
        non_zero_params={'clu':val}
    elif par == 72:
        non_zero_params={'cld':val}
    elif par == 73:
        non_zero_params={'cqe':val}
    elif par == 74:
        non_zero_params={'cqu1':val}
    elif par == 75:
        non_zero_params={'cqu8':val}
    elif par == 76:
        non_zero_params={'cqd1':val}
    elif par == 77:
        non_zero_params={'cqd8':val}
    elif par == 78:
        non_zero_params={'cledqAbs':val}
    elif par == 79:
        non_zero_params={'cquqd1Abs':val}
    elif par == 80:
        non_zero_params={'cquqd8Abs':val}
    elif par == 81:
        non_zero_params={'clequ1Abs':val}
    elif par == 82:
        non_zero_params={'clequ3Abs':val}

    else: 
        raise RuntimeError("runNumber %i not recognised in these jobOptions."%par)

    return non_zero_params
