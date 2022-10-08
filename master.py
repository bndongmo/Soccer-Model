# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 11:04:31 2020

@author: bryan

"""
e = 2.718281828459045
up = 1.5
down = 1.5**(1/3) 
#IMPORTS
import pandas 
import numpy
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
import time
import sympy as sym
from sympy import nsolve
import math
import shelve
import profile

#FUNCTIONS
memo = {}
def factorial(x):
    if x in memo:
        return memo[x]
    if x == 0 or x == 1:
        f = 1
    else:
        f = x * factorial(x - 1)
    memo[x] = f
    return f

comb_memo = {}
def comb(x, y):
    if y == 0 or x == y:
        return 1
    if y == 1 or y == x-1:
        return x
    if (x, y) in comb_memo:
        return comb_memo[(x,y)]
    c = factorial(x)/(factorial(y)*factorial(x-y))
    comb_memo[(x, y)] = c
    return c

mult_memo = {}
def multiply(x, y):
    if x == 1:
        return y
    if y == 1:
        return x
    if x == 0 or y == 0:
        return 0
    
    if (x, y) in mult_memo:
        return mult_memo[(x, y)]
    m = x * y
    mult_memo[(x, y)] = m
    return m

exp_memo = {}
def exp(x, y):
    if y == 0:
        return 1
    if y == 1:
        return x
    if (x, y) in exp_memo:
        return exp_memo[(x, y)]
    ex = x ** y
    exp_memo[(x, y)] = ex
    return ex

div_memo = {}
def divide(x, y):
    if y == 1:
        return x
    if x == y:
        return 1
    if (x, y) in div_memo:
        return div_memo[(x, y)]
    d = x / y
    div_memo[(x, y)] = d
    
    return d


prob_memo = {}

def probs(z, home_mu, away_mu, p, sh=False):
    with shelve.open('probs') as prob_memo_shelf:
        if str((z, home_mu, away_mu, p)) in prob_memo_shelf:
            return prob_memo_shelf[str((z, home_mu, away_mu, p))]
        
        if str((z, home_mu, away_mu, p)) in prob_memo:
            return prob_memo[str((z, home_mu, away_mu, p))]
        probs = {}
        o = 0
        u = 0
        total = 0                
        for x in range(10): 
            for y in range(10):
                prob = 0
                
                coefficient = divide(divide(multiply(multiply(exp(e, -home_mu-away_mu-z), exp(home_mu, x)), exp(away_mu, y)), factorial(x)), factorial(y))    
                for i in range(min(x, y)+1):
                    prob += multiply(multiply(multiply(comb(x, i), comb(y, i)), factorial(i)), exp(divide(z, multiply(home_mu, away_mu)), i))    
                prob = multiply(prob, coefficient)
                #prob = multiply(prob, 1-p)
        
                probs[(x, y)] = prob
                total += prob
                
               
        
        
        o_draw = 0
        for key in probs:
            probs[key] = multiply(probs[key], divide(1, total))
            probs[key] = multiply(probs[key], 1-p)
            if key[0] == key[1]:
                o_draw += probs[key]
                
        ahome = 0 #probabilities
        aaway = 0
        adraw = 0
        
        for key in probs:
            x = key[0]
            y = key[1]
            prob = probs[key]
            if x > y:
                ahome += prob
                if x+y > 2.5:
                    o += prob
                else:
                    u += prob
            elif x < y:
                aaway += prob
                if x+y > 2.5:
                    o += prob
                else:
                    u += prob
            else:
                adjust = multiply(divide(prob, o_draw), p)
                adraw += prob + adjust
                probs[key] += adjust
                if x+y > 2.5:
                    o += prob + adjust
                else:
                    u += prob + adjust
                    
        prob_memo[str((z, home_mu, away_mu, p))] = probs
        if sh:
            prob_memo_shelf[str((z, home_mu, away_mu, p))] = probs
        return probs

def binary_search(arr, val, start, end):
     
    # we need to distinugish whether we
    # should insert before or after the
    # left boundary. imagine [0] is the last
    # step of the binary search and we need
    # to decide where to insert -1
    if start == end:
        if arr[start][1][3] > val:
            return start
        else:
            return start+1
 
    # this occurs if we are moving
    # beyond left's boundary meaning
    # the left boundary is the least
    # position to find a number greater than val
    if start > end:
        return start
 
    mid = (start+end)//2
    if arr[mid][1][3] < val:
        return binary_search(arr, val, mid+1, end)
    elif arr[mid][1][3] > val:
        return binary_search(arr, val, start, mid-1)
    else:
        return mid
    
def extract_hcap(fav, typ, hcap, output):   
    if not typ:
        if abs(hcap).is_integer():
            typ = 'whole'
        else:
            typ = 'half'
      
    if hcap == 0:
        win = output[0]
        push = output[1]
        lose = output[2]
        return (multiply(win, divide(1, 1-push)), multiply(lose, divide(1, 1-push)))
        
    if typ == 'whole':
        if fav == 'home':
            if hcap == -1:
                win = output[5]
                push = output[0]-output[5]
                lose = 1 - win - push                                               
            elif hcap == -2:
                win = output[6]
                push = output[5]-output[6]
                lose = 1 - win - push               
            elif hcap == -3:
                win = output[7]
                push = output[6]-output[7]
                lose = 1 - win - push                     
        else:
            if hcap == 1:
                lose = output[8]
                push = output[2]-output[8]
                win = 1 - lose - push               
            elif hcap == 2:
                lose = output[9]
                push = output[8]-output[9]
                win = 1 - lose - push                
            elif hcap == 3:
                lose = output[10]
                push = output[9]-output[10]
                win = 1 - lose - push               
        return (multiply(win, divide(1, 1-push)), multiply(lose, divide(1, 1-push)))
    
    if typ == 'half':
        if fav == 'home':
            if hcap == -0.5:
                win = output[0]
                lose = 1 - win 
            elif hcap == -1.5:
                win = output[5]
                lose = 1 - win               
            elif hcap == -2.5:
                win = output[6]
                lose = 1 - win
                
            elif hcap == -3.5:
                win = output[7]
                lose = 1 - win
        else:    
            if hcap == 0.5:
                lose = output[2]
                win = 1 - lose
            elif hcap == 1.5:
                lose = output[8]
                win = 1 - lose
                
            elif hcap == 2.5:
                lose = output[9]
                win = 1 - lose
                
            elif hcap == 3.5:
                lose = output[10]
                win = 1 - lose      

        return (win, lose)
    
    else:
        if abs(hcap - .25).is_integer():
            a = extract_hcap(fav, None, hcap-.25, output)
            b = extract_hcap(fav, None, hcap+.25, output)
        else:
            a = extract_hcap(fav, None, hcap-.25, output)
            b = extract_hcap(fav, None, hcap+.25, output)
            
        return (divide(a[0]+b[0], 2), divide(a[1]+b[1], 2))
    
def parameters(home_prob, draw_prob, away_prob, o_prob, u_prob, hcap, pahhodds, pahaodds):
    min_score = float('inf')
    sol = None
    start = binary_search(split, o_prob-.02, 0, len(split)-1)
    end = binary_search(split, o_prob+.02, 0, len(split)-1)
    
    for i in range(start, end):
        #if split[i][0][0] == 0 and split[i][0][3] == 0: 
        score = 0
        score += divide(abs(home_prob-split[i][1][0]), 3)
        score += divide(abs(draw_prob-split[i][1][1]), 3)
        score += divide(abs(away_prob-split[i][1][2]), 3)
        score += divide(abs(o_prob-split[i][1][3]), 2)
        score += divide(abs(u_prob-split[i][1][4]), 2)
        if score < min_score:
            min_score = score
            sol = split[i]
            
    min_score = float('inf')
          
    if hcap < 0:
        fav = 'home'
    elif hcap > 0:
        fav = 'away'
    else:
        fav = None
    if abs(hcap).is_integer():
        typ = 'whole'
    elif abs(hcap*2).is_integer():
        typ = 'half'
    else:
        typ = 'quarter'
    _z = sol[0][0]
    _hmu = sol[0][1]
    _amu = sol[0][2]
    _p = sol[0][3]
    for a in range(int(multiply(max(0, _z-.05), 100)), int(multiply(_z+.05,100)), 2):
        for b in range(int(multiply(_hmu-.05,100)), int(multiply(_hmu+.05,100)), 2):
            for c in range(int(multiply(_amu-.05,100)), int(multiply(_amu+.05,100)), 2):
                for d in range(int(multiply(max(0, _p-.005), 1000)), int(multiply(_p+.005,1000)), 2):
                    predict= [[divide(a,100), divide(b,100), divide(c,100), divide(d,1000)]]
                    predict_ = poly.fit_transform(predict)
                    candidate = clf.predict(predict_)
                    ahprobs = extract_hcap(fav, typ, hcap, candidate[0])
                    score = 0
                    score += divide(abs(home_prob-candidate[0][0]), 6)
                    score += divide(abs(draw_prob-candidate[0][1]), 6)
                    score += divide(abs(away_prob-candidate[0][2]), 6)
                    score += divide(abs(o_prob-candidate[0][3]), 2)
                    score += divide(abs(u_prob-candidate[0][4]), 2)
                    score += divide(abs(pahhodds-ahprobs[0]), 4)
                    score += divide(abs(pahaodds-ahprobs[1]), 4)
                    if score < min_score:
                        min_score = score
                        sol = ((divide(a,100), divide(b, 100), divide(c, 100), divide(d, 1000)), (candidate[0][0], candidate[0][1], candidate[0][2], candidate[0][3], candidate[0][4]))
                        
    _z = sol[0][0]
    _hmu = sol[0][1]
    _amu = sol[0][2]
    _p = sol[0][3]
       
    for a in range(int(multiply(max(0, _z-.02), 100)), int(multiply(_z+.02,100))):
        for b in range(int(multiply(_hmu-.02,100)), int(multiply(_hmu+.02,100))):
            for c in range(int(multiply(_amu-.02,100)), int(multiply(_amu+.02,100))):
                for d in range(int(multiply(max(0, _p-.002), 1000)), int(multiply(_p+.002,1000))):
                    predict= [[divide(a,100), divide(b,100), divide(c,100), divide(d,1000)]]
                    predict_ = poly.fit_transform(predict)
                    candidate = clf.predict(predict_)
                    ahprobs = extract_hcap(fav, typ, hcap, candidate[0])
                    score = 0
                    score += divide(abs(home_prob-candidate[0][0]), 6)
                    score += divide(abs(draw_prob-candidate[0][1]), 6)
                    score += divide(abs(away_prob-candidate[0][2]), 6)
                    score += divide(abs(o_prob-candidate[0][3]), 2)
                    score += divide(abs(u_prob-candidate[0][4]), 2)
                    score += divide(abs(pahhodds-ahprobs[0]), 4)
                    score += divide(abs(pahaodds-ahprobs[1]), 4)
                    if score < min_score:
                        min_score = score
                        sol = ((divide(a,100), divide(b, 100), divide(c, 100), divide(d, 1000)), (candidate[0][0], candidate[0][1], candidate[0][2], candidate[0][3], candidate[0][4]))
    return sol[0]

def func(x):
    return (x**2 + 2*x + 2) / (2*e**x)

memo_2 = {}
def find_mu(o, u):
    if (o, u) in memo_2:
        return memo_2[(o, u)]
    overround = 1/o + 1/u
    u_prob = 1/u/overround
    epsilon = .001
    mu = 1.5
    while True:
        if abs(func(mu) - u_prob) <= epsilon:
            memo_2[(o, u)] = mu
            return mu
        mu += .001
        
def goal_lines(prec_line):
    if prec_line >= 4.162:
        rounded_line = 4.5
    elif prec_line >= 3.16:
        rounded_line = 3.5
    elif prec_line >= 2.156:
        rounded_line = 2.5
    elif prec_line >= 1.146:
        rounded_line = 1.5
    else:
        rounded_line = 0.5
        
    if prec_line >= 4.522:
        rounded_line_2 = 4.5
    elif prec_line >= 4.272:
        rounded_line_2 = 4.25
    elif prec_line >= 4.053:
        rounded_line_2 = 4
    elif prec_line >= 3.813:
        rounded_line_2 = 3.75
    elif prec_line >= 3.518:
        rounded_line_2 = 3.5
    elif prec_line >= 3.268:
        rounded_line_2 = 3.25
    elif prec_line >= 3.053:
        rounded_line_2 = 3
    elif prec_line >= 2.816:
        rounded_line_2 = 2.75
    elif prec_line >= 2.512:
        rounded_line_2 = 2.5
    elif prec_line >= 2.262:
        rounded_line_2 = 2.25
    elif prec_line >= 2.052:
        rounded_line_2 = 2
    elif prec_line >= 1.822:
        rounded_line_2 = 1.75
    else:
        rounded_line_2 = 1.5
        
    return float(rounded_line), float(rounded_line_2)

def prob_breakdown(probs, home_line, away_line, goal_line, hcap, goal_line_2, hcap_2):
    under_prob = 0
    under_prob_2 = 0
    under_push_prob = 0
    under_half_win_prob = 0
    under_half_lose_prob = 0
    home_under_prob = 0
    away_under_prob = 0
    home_under_push_prob = 0
    away_under_push_prob = 0
    home_cover_prob = 0
    home_cover_prob_2 = 0
    home_push_prob = 0
    home_half_win_prob = 0
    home_half_lose_prob = 0
    home_win_prob = 0
    draw_prob = 0
    btts_yes = 0
    
    home_and_btts_yes = 0
    home_to_nil = 0
    home_cs = 0
    home_or_btts = 0
    
    draw_and_btts_yes = 0
    draw_and_btts_no = 0
    draw_or_btts = 0
    
    away_and_btts_yes = 0
    away_to_nil = 0
    away_cs = 0
    away_or_btts = 0
    
    for key in probs:
        i = key[0]
        j = key[1]
        occurence_prob = probs[key]
        #game total
        if i + j <= goal_line - 0.5:
            under_prob += occurence_prob
            
        #game total 2
        if i + j <= goal_line_2 - 0.5:
            under_prob_2 += occurence_prob
        elif i + j - goal_line_2 == .25:
            under_half_lose_prob += occurence_prob
        elif i + j - goal_line_2 == -.25:
            under_half_win_prob += occurence_prob
        elif i + j == goal_line_2:   
            under_push_prob += occurence_prob 
            
        #home total
        if i < home_line:
            home_under_prob += occurence_prob
        elif i == home_line:
            home_under_push_prob += occurence_prob
            
        #away total
        if j <= away_line - 0.5:
            away_under_prob += occurence_prob
        elif j == away_line:
            away_under_push_prob += occurence_prob
            
        #hcap
        if j - i <= hcap - 0.5:
            home_cover_prob += occurence_prob
       
        #hcap2
        if j - i <= hcap_2 - 0.5:
            home_cover_prob_2 += occurence_prob
        elif i + hcap_2 == j - .25:
            home_half_lose_prob += occurence_prob
        elif i + hcap_2 == j + .25:
            home_half_win_prob += occurence_prob
        elif j - i == hcap_2:
            home_push_prob += occurence_prob
        
        home = False
        draw = False
        away = False
        #result
        if i > j:
            home_win_prob += occurence_prob
            home = True
        elif i == j:
            draw_prob += occurence_prob
            draw = True
        else:
            away = True
        
        btts = False
        #btts
        if i > 0 and j > 0:
            btts_yes += occurence_prob
            btts = True
            
        #resultand
        if home:
            if btts:
                home_and_btts_yes += occurence_prob
            else:
                home_to_nil += occurence_prob
        elif draw:
            if btts:
                draw_and_btts_yes += occurence_prob
            else:
                draw_and_btts_no += occurence_prob
        elif away:
            if btts:
                away_and_btts_yes += occurence_prob
            else:
                away_to_nil += occurence_prob
                
        #resultor
        if home or btts:
            home_or_btts += occurence_prob
        if draw or btts:
            draw_or_btts += occurence_prob
        if away or btts:
            away_or_btts += occurence_prob
            
        #cs
        if j == 0:
            home_cs += occurence_prob
        if i == 0:
            away_cs += occurence_prob
            
                
    over_prob = 1 - under_prob
    over_prob_2 = 1 - under_prob_2 - under_push_prob - under_half_win_prob - under_half_lose_prob
    home_over_prob = 1 - home_under_prob - home_under_push_prob
    away_over_prob = 1 - away_under_prob - away_under_push_prob
    away_cover_prob = 1 - home_cover_prob
    away_cover_prob_2 = 1 - home_cover_prob_2 - home_push_prob - home_half_win_prob - home_half_lose_prob
    away_win_prob = 1 - home_win_prob - draw_prob
    btts_no = 1 - btts_yes
    
    return {'result': [home_win_prob, draw_prob, away_win_prob],
            'hcap': [home_cover_prob, 0, away_cover_prob],
            'hcap_2': [home_cover_prob_2, home_half_win_prob, home_push_prob, home_half_lose_prob, away_cover_prob_2],
            'game_total': [over_prob, 0, under_prob],
            'game_total_2': [over_prob_2, under_half_lose_prob, under_push_prob, under_half_win_prob, under_prob_2],
            'home_total': [home_over_prob, home_under_push_prob, home_under_prob],
            'away_total': [away_over_prob, away_under_push_prob, away_under_prob],
            'btts': [btts_yes, btts_no],
            'win+btts': [home_and_btts_yes, away_and_btts_yes],
            'draw+btts': [draw_and_btts_yes, draw_and_btts_no],
            'home_to_nil': [home_to_nil, 1-home_to_nil],
            'away_to_nil': [away_to_nil, 1-away_to_nil],
            'home_cs': [home_cs, 1 - home_cs],
            'away_cs': [away_cs, 1 - away_cs],
            'home_or_btts': [home_or_btts, 1 - home_or_btts],
            'draw_or_btts': [draw_or_btts, 1 - draw_or_btts],
            'away_or_btts': [away_or_btts, 1 - away_or_btts]}

def calc_kelly(odds, win_prob, push_prob = 0):
    lose_prob = 1 - win_prob - push_prob
    return (win_prob * (odds-1) - lose_prob)/((lose_prob + win_prob)*(odds-1))

def nearest(hcap, h, a):
    l = [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5]
    for line in l:
        if abs(line-hcap) <= .25:
            return line
        if abs(line-hcap) == .5:
            if h < a:           
                return line
            return line+1
        
def kelly_hw(odds, w, hw):
    if hw == 0:
        return calc_kelly(odds, w)
    
    m = odds-1
    num1 = m*hw+2*hw-3+m*w+w
    num2 = (m**2*hw**2+4*hw**2+4*m*hw**2+4*w*hw+2*m**2*w*hw+6*m*w*hw-4*hw-2*m*hw+m**2*w**2+2*m*w**2+2*m*w+w**2+2*w+1)**(1/2)
    den = 2 * m

    return (num1+num2)/den
    

def kelly_hl(odds, w, hl):
    if hl == 0:
       return calc_kelly(odds, w)
    m = odds-1
    num1 = -(m*hl+1-m*w-w-2*m)
    num2 = (m**2*w**2+2*m*w**2+w**2-4*m**2*w-6*m*w-2*m**2*w*hl-2*m*w*hl-2*w+4*m**2+4*m+m**2*hl**2-2*m*hl-4*m**2*hl+1)**(1/2)
    den = 2 * m
    return (num1-num2)/den

    
    

def trim_list(list_):
    if len(list_) > 8:
        list_.pop(0)
        
def weighted_total(l):
    n = 0
    for i in range(len(l)):
        n += l[i] * (len(l) - i) 
    return n / 36 * 8

def insert(pick, picks):
    if not picks:
        picks.append(pick)
        
    else:
        score = pick[0][1]
        low = 0
        high = len(picks)-1
        while True: 
            mid = (low + high) // 2
            if (picks[mid][0][1] == score):
                index = mid
                break
            
            elif (picks[mid][0][1] < score): 
              low = mid + 1
              if (low > high):
                index = mid + 1
                break
            else: 
                high = mid - 1
                if (low > high):
                  index = mid
                  break
        picks = picks[0:index] + [pick] + picks[index:]
            
    return picks

def res(pick, home_g, away_g):
    plus = pick[3] >= 2
    if plus:
        win = pick[3] - 1
        lose = -1
    else:
        win = 1
        lose = -(1/(pick[3]-1))
        
    if home in pick[0]:
        pro = home_g
        con = away_g
    else:
        pro = away_g
        con = home_g
    double = 'vs' in pick[0] 
    #3WAY
    if "ML" in pick[0]:
        if pro > con:
            return win
        return lose
    if "Draw " in pick[0]:
        if pro == con:
            return win
        return lose
    #TOTALS
    if "O-" in pick[0]:
        if double:
            if home_g+away_g >= pick[2] + .5:
                return win
            if home_g+away_g == pick[2] + .25:
                return win/2
            if home_g+away_g == pick[2]:
                return 0
            if home_g+away_g == pick[2] - .25:
                return lose/2
            return lose
        else:
            if pro > pick[2]:
                return win
            return lose
    if "U-" in pick[0]:
        if double:
            if home_g+away_g >= pick[2] + .5:
                return lose
            if home_g+away_g == pick[2] + .25:
                return lose/2
            if home_g+away_g == pick[2]:
                return 0
            if home_g+away_g == pick[2] - .25:
                return win/2
            return win
        else:
            if pro > pick[2]:
                return lose
            return win
        
    #BTTS/0.5/RESULT+BTTS
    if "BTTS" in pick[0]:
        if "DOBTTSY" in pick[0]:
            if pro == con or (pro > .5 and con >.5):
                return win
            return lose
        if "DOBTTSN" in pick[0]:
            if pro != con and (pro < .5 or con < .5):
                return win
            return lose
        
        if "OBTTSY-" in pick[0]:
            if pro > 0.5:
                return win
            return lose
        if "OBTTSN-" in pick[0]:
            if pro < 0.5:
                return win
            return lose     
        
        if "HBTTSY" in pick[0] or "ABTTSY" in pick[0]:
            if pro > con and con > .5:
                return win
            return lose
        
        if "DBTTSY" in pick[0]:
            if pro == con and con > .5:
                return win
            return lose
        
        if "DBTTSN" in pick[0]:
            if pro == 0 and con == 0:
                return win
            return lose
        
        if "Yes" in pick[0]:
            if pro > .5 and con > .5:
                return win
            return lose
        
        if pro < .5 or con < .5:
            return win
        return lose
    
    if "HTNY" in pick[0]:
        if home_g > away_g and away_g < .5:
            return win
        return lose
    
    if "ATNY" in pick[0]:
        if home_g < away_g and home_g < .5:
            return win
        return lose
    
    if "HTNN" in pick[0]:
        if home_g <= away_g or away_g > .5:
            return win
        return lose
    
    if "ATNN" in pick[0]:
        if home_g >= away_g or home_g > .5:
            return win
        return lose
    
    #DNB/ANB/HNB/DC
    if "DNB-" in pick[0]:
        if pro > con:
            return win
        if pro == con:
            return 0
        return lose
    if "DC-" in pick[0]:
        if pro >= con:
            return win
        return lose
    if "ANB-" in pick[0]:
        if "DRAW" in pick[0]:
            if pro == con:
                return win
            if away_g > home_g:
                return 0
            return lose
        if pro > con:
            return win
        if away_g > home_g:
            return 0
        return lose
    
    if "HNB-" in pick[0]:
        if "DRAW" in pick[0]:
            if pro == con:
                return win
            if away_g < home_g:
                return 0
            return lose
        if pro > con:
            return win
        if away_g < home_g:
            return 0
        return lose
    actual = pro - con + pick[2]
    if actual >= .5:
        return win
    elif actual == .25:
        return win / 2
    elif actual == 0:
        return 0
    elif actual == -.25:
        return lose / 2
    else:
        return lose
    
def maxSubArray(arr):
    if not arr:
        return 0
    if len(arr) == 1:
        return arr[0]
    currentSum = arr[0]     #Starting with the first element
    maxSum = currentSum
    for i in range(1, len(arr)):
      currentSum = max(arr[i], currentSum + arr[i])     #Updating current sum if the current element adds more sum value
      maxSum = max(maxSum, currentSum)                  #Updating maxSum
    return maxSum    

def plus(odds):
    return 1 + odds/100

def minus(odds):
    return 1 + 100/odds 

def sort(picks, i = 0, center = 1000):
    if i < len(picks):
        max_ = -float('inf')
        index = None
        for j in range(i, len(picks)):
            if 100 - abs(picks[j][1]-center) > max_:
                max_ = 100 - abs(picks[j][1]-center)
                index = j
        temp = picks[i]
        picks[i] = picks[index]
        picks[index] = temp
        sort(picks, i+1, center)
           
#FILES
with open("results.txt", encoding="utf-8") as f:
        results = f.readlines()
        
split = []
for r in results:
    split.append(r.split("     "))
 
for s in split:
    s[0] = eval(s[0])
    s[1] = eval(s[1])
    
split.reverse()

X = []
vector = []
with open("regress.txt", encoding="utf-8") as f:     
        for line in f:
            line = f.readline().split('\t')
            X.append(line[:4])
            for i in range(len(X[-1])):
                X[-1][i] = float(X[-1][i])
            vector.append(line[4:])
            for i in range(len(vector[-1])):
                vector[-1][i] = float(vector[-1][i])
            


 
poly = PolynomialFeatures(degree=10)
X_ = poly.fit_transform(X)
clf = linear_model.LinearRegression()
clf.fit(X_, vector) 
        
leagues = []
fixtures = []
with open("family.txt", encoding="utf-8") as f:
        leagues.append(f.readlines())
with open("fixtures.txt", encoding="utf-8") as f:
        fixtures = f.readlines()
m = []
n = []
for league in leagues:
    for line in league:
        result = line.split("\t")
        m.append(result)

for fixture in fixtures:
    n.append(fixture.split('\t'))
    

fixtures = n


team_scores = {}
defense_scores = {}
xgs = {}
xgs_ag = {}
lines = {}
lines_ag = {}
team_lines = {}   
opp_team_lines = {}   
total_goals = {}
total_goals_ag = {}
total_xg = {}
total_xg_ag = {}
with shelve.open('pinn_params') as pinn_params_memo:
    with shelve.open('avg_params') as avg_params_memo:
        for game in m:
            if game[3] != 'HomeTeam': 
                home = game[3]
                away = game[4]
                date = game[1]
                try:
                    home_cover_odds = float(game[61])
                    away_cover_odds = float(game[62])
                except:
                    try:
                        home_cover_odds = float(game[59])
                        away_cover_odds = float(game[60])
                    except:
                         home_cover_odds = float(game[65])
                         away_cover_odds = float(game[66])
                        
                
                try:
                    o = float(game[52])
                    u = float(game[53])
                except:
                    try:
                        o = float(game[50])
                        u = float(game[51])
                    except:
                        o = float(game[56])
                        u = float(game[57])
                try:
                    home_odds = float(game[35])
                    draw_odds = float(game[36])
                    away_odds = float(game[37])
                except:
                    try:
                        home_odds = float(game[26])
                        draw_odds = float(game[27])
                        away_odds = float(game[28])
                    except:
                        home_odds = float(game[47])
                        draw_odds = float(game[48])
                        away_odds = float(game[49])
                home_g = float(game[5])
                away_g = float(game[6])
                home_xg = float(game[7])
                away_xg = float(game[8])
                hcap = float(game[-9])
                line = find_mu(o, u)
                overround = 1/home_odds+1/draw_odds+1/away_odds
                home_prob = 1/home_odds/overround
                draw_prob = 1/draw_odds/overround
                away_prob = 1/away_odds/overround
                overround = 1/o+1/u
                o_prob = 1/o/overround
                u_prob = 1/u/overround
                overround = 1/home_cover_odds + 1/away_cover_odds
                pahhodds = 1/home_cover_odds/overround
                pahaodds = 1/away_cover_odds/overround 
                if str((home, away, date)) in pinn_params_memo:
                    pinn_params = pinn_params_memo[str((home, away, date))]
                else:
                    pinn_params = parameters(home_prob, draw_prob, away_prob, o_prob, u_prob, hcap, pahhodds, pahaodds)    
                    pinn_params_memo[str((home, away, date))] = pinn_params
        
                z, home_mu, away_mu, p = pinn_params
                home_line = home_mu+z/2
                away_line = away_mu+z/2
                lines[game[3]] = lines.setdefault(game[3], []) + [home_line]
                trim_list(lines[game[3]])
                lines[game[4]] = lines.setdefault(game[4], []) + [away_line]
                trim_list(lines[game[4]])
                lines_ag[game[3]] = lines_ag.setdefault(game[3], []) + [away_line]
                trim_list(lines_ag[game[3]])
                lines_ag[game[4]] = lines_ag.setdefault(game[4], []) + [home_line]
                trim_list(lines_ag[game[4]])
                xgs[game[3]] = xgs.setdefault(game[3], []) + [home_xg]
                trim_list(xgs[game[3]])
                xgs[game[4]] = xgs.setdefault(game[4], []) + [away_xg]
                trim_list(xgs[game[4]])
                xgs_ag[game[3]] = xgs_ag.setdefault(game[3], []) + [away_xg]
                trim_list(xgs_ag[game[3]])
                xgs_ag[game[4]] = xgs_ag.setdefault(game[4], []) + [home_xg]
                trim_list(xgs_ag[game[4]])
                
                home_goal_score = home_xg/home_line
                away_goal_score = away_xg/away_line
                if home_goal_score > 10:
                    print('dsfsf')
                    home_goal_score = 1
                if away_goal_score > 10:
                    print('sfssfsdf')
                    away_goal_score = 1
                away_def_score = home_goal_score
                home_def_score = away_goal_score
        
                team_scores[game[3]] = team_scores.setdefault(game[3], []) + [home_goal_score]
                trim_list(team_scores[game[3]])
                team_scores[game[4]] = team_scores.setdefault(game[4], []) + [away_goal_score]
                trim_list(team_scores[game[4]])
                defense_scores[game[3]] = defense_scores.setdefault(game[3], []) + [home_def_score]
                trim_list(defense_scores[game[3]])
                defense_scores[game[4]] = defense_scores.setdefault(game[4], []) + [away_def_score]
                trim_list(defense_scores[game[4]])
                
           
                total_goals[game[3]] = total_goals.setdefault(game[3], 0) + home_g
                total_goals[game[4]] = total_goals.setdefault(game[4], 0) + away_g
                total_goals_ag[game[4]] = total_goals_ag.setdefault(game[4], 0) + home_g
                total_goals_ag[game[3]] = total_goals_ag.setdefault(game[3], 0) + away_g
                
                total_xg[game[3]] = total_xg.setdefault(game[3], 0) + home_xg
                total_xg[game[4]] = total_xg.setdefault(game[4], 0) + away_xg
                total_xg_ag[game[4]] = total_xg_ag.setdefault(game[4], 0) + home_xg
                total_xg_ag[game[3]] = total_xg_ag.setdefault(game[3], 0) + away_xg
        
lethalness = {}
gk = {}
for key in total_xg:
    lethalness[key] = total_goals[key] / total_xg[key]  
    gk[key] = total_goals_ag[key] / total_xg_ag[key]  

                    
            

#xgs['Moreirense'] = []
ratings = {}
goal_ratings = {}
team_goal_ratings = {}
def_ratings = {}
goal_stdevs = {}

for key in team_scores.keys():
    team_goal_ratings[key] = weighted_total(xgs[key]) / weighted_total(lines[key])#weighted_geomean(team_scores[key])#sum(team_scores[key]) / len(team_scores[key])
    def_ratings[key] = weighted_total(xgs_ag[key]) / weighted_total(lines_ag[key])#weighted_geomean(defense_scores[key])#sum(defense_scores[key]) / len(defense_scores[key])

#print(team_goal_ratings)
delete = [key for key in team_goal_ratings if len(team_scores[key]) < 8] 
# delete the key 
for key in delete: del team_goal_ratings[key] 
for key in delete: del def_ratings[key] 
#print(team_goal_ratings['Newcastle'])
print(max(team_goal_ratings, key=lambda key: team_goal_ratings[key]))
print(min(team_goal_ratings, key=lambda key: team_goal_ratings[key]))
print(max(def_ratings, key=lambda key: def_ratings[key]))
print(min(def_ratings, key=lambda key: def_ratings[key]))    
blacklist = []
check = []



betlist = []
w = 2
x = 4
y = 4
z = 1
letters = {0: (1, 2, 1, 0), 1: (1, 2, 1, 2), 2: (3, 0, 1, 0)}
picks = []
dog_picks = []
fav_picks = []
one = []
two = []
betlists = {}
aways = []
with shelve.open('pinn_params') as pinn_params_memo:
    with shelve.open('avg_params') as avg_params_memo:
        for let in range(3):
            w = letters[let][0]
            x = letters[let][1]
            y = letters[let][2]
            z = letters[let][3]
            for game in fixtures:
                #print(fixture)
                #os.system('clear')
                if fixture[3]:#== '4/20/2021' :#and (fixture[0] == 'E0' or fixture[0] == 'D1' or fixture[0] == 'SP1' or fixture[0] == 'I1' or fixture[0] == 'F1'):
                    if len(team_scores[game[3]]) >= 8 and len(team_scores[game[4]]) >= 8:
                        home = game[3]
                        away = game[4]    
                        date = game[1]
                        hcap = float(game[43])
                        orig = hcap
                        home_odds = float(game[32])
                        draw_odds = float(game[33])
                        away_odds = float(game[34])
                        home_cover_odds = float(game[50]) 
                        away_cover_odds = float(game[51]) 
                        orig_home_cover_odds = home_cover_odds
                        orig_away_cover_odds = away_cover_odds
                        o = float(game[41])
                        u = float(game[42])
                        
                        try:
                            pinn_home_cover_odds = float(game[46])
                            pinn_away_cover_odds = float(game[47])
                        except:
                            try:
                                pinn_home_cover_odds = float(game[44])
                                pinn_away_cover_odds = float(game[45])
                            except:
                                pinn_home_cover_odds = home_cover_odds
                                pinn_away_cover_odds = away_cover_odds
                                
                        
                        try:
                            pinn_o = float(game[37])
                            pinn_u = float(game[38])
                        except:
                            try:
                                pinn_o = float(game[35])
                                pinn_u = float(game[36])
                            except:
                                pinn_o = o
                                pinn_u = u
                        try:
                            pinn_home_odds = float(game[20])
                            pinn_draw_odds = float(game[21])
                            pinn_away_odds = float(game[22])
                        except:
                            try:
                                pinn_home_odds = float(game[11])
                                pinn_draw_odds = float(game[12])
                                pinn_away_odds = float(game[13])
                            except:
                                pinn_home_odds = home_odds
                                pinn_draw_odds = draw_odds
                                pinn_away_odds = away_odds
                        
                        prec_line = find_mu(o, u)  
                        
                        overround = 1/pinn_home_odds+1/pinn_draw_odds+1/pinn_away_odds
                        home_prob = 1/pinn_home_odds/overround
                        draw_prob = 1/pinn_draw_odds/overround
                        away_prob = 1/pinn_away_odds/overround
                        overround = 1/pinn_o+1/pinn_u
                        o_prob = 1/pinn_o/overround
                        u_prob = 1/pinn_u/overround
                        overround = 1/pinn_home_cover_odds + 1/pinn_away_cover_odds
                        pahhodds = 1/pinn_home_cover_odds/overround
                        pahaodds = 1/pinn_away_cover_odds/overround 
                        
                        if str((home, away, date)) in pinn_params_memo:
                            pinn_params = pinn_params_memo[str((home, away, date))]
                        else:
                            pinn_params = parameters(home_prob, draw_prob, away_prob, o_prob, u_prob, hcap, pahhodds, pahaodds)    
                            pinn_params_memo[str((home, away, date))] = pinn_params
            
                        z, home_mu, away_mu, p = pinn_params
                        home_line = home_mu+z/2
                        away_line = away_mu+z/2
                        
                        bets = []
                        overround = 1/home_odds+1/draw_odds+1/away_odds
                        home_prob = 1/home_odds/overround
                        draw_prob = 1/draw_odds/overround
                        away_prob = 1/away_odds/overround
                        overround = 1/o+1/u
                        o_prob = 1/o/overround
                        u_prob = 1/u/overround
                        overround = 1/home_cover_odds + 1/away_cover_odds
                        pahhodds = 1/home_cover_odds/overround
                        pahaodds = 1/away_cover_odds/overround 
                        hcap = nearest(hcap, home_cover_odds, away_cover_odds)
                        if str((home, away, date)) in avg_params_memo:
                            avg_params = avg_params_memo[str((home, away, date))]
                        else:
                            avg_params = parameters(home_prob, draw_prob, away_prob, o_prob, u_prob, hcap, pahhodds, pahaodds)    
                            avg_params_memo[str((home, away, date))] = avg_params
                  
                        rounded_line, rounded_line_2 = goal_lines(prec_line)
                        #rounded_line adj
                        
                      
                        home_mult = max((team_goal_ratings[home])**(w/4) * (def_ratings[away])**(x/4) * (lethalness[home])**(y/4) * (gk[away])**(z/4), 0.01)
                        away_mult = max((team_goal_ratings[away])**(w/4) * (def_ratings[home])**(x/4) * (lethalness[away])**(y/4) * (gk[home])**(z/4), 0.01)
                    
                    
                        
                        z, home_mu, away_mu, p = avg_params
                        avg_home_line = home_mu+z/2
                        avg_away_line = away_mu+z/2
                        hrounded_line, arounded_line = goal_lines(avg_home_line)[0], goal_lines(avg_away_line)[0]
                        #line adj
                    
                            
 
                        other_bdowns = prob_breakdown(probs(z, home_mu, away_mu, p, True), hrounded_line, arounded_line, rounded_line, hcap, rounded_line_2, orig)
                        z, home_mu, away_mu, p = pinn_params
                        scores = probs(z, home_line*home_mult-z/2, away_line*away_mult-z/2, p, True)
                        breakdowns = prob_breakdown(scores, hrounded_line, arounded_line, rounded_line, hcap, rounded_line_2, orig) 
                     
                        
                        if home_mult >= .1 and away_mult >= .1: 
                            #3WAY
                    
                            home_win = breakdowns['result'][0]
                            draw = breakdowns['result'][1]
                            away_win = breakdowns['result'][2]

                            kelly = calc_kelly(home_odds, home_win)
                            bets.append((home + " ML", kelly, orig, home_odds))
           
                            kelly = calc_kelly(away_odds, away_win)
                            bets.append((away + " ML", kelly, -orig, away_odds))
        
                            kelly = calc_kelly(draw_odds, draw)
                            bets.append(("Draw " + home + " vs " + away, kelly, orig, draw_odds))
                            
                            #HANDICAPS
                            home_cover_odds = 1 / (other_bdowns['hcap'][0] / (1 - other_bdowns['hcap'][1]) * 1.07)
                            away_cover_odds = 1 / (other_bdowns['hcap'][2] / (1 - other_bdowns['hcap'][1]) * 1.07)
                            home_cover = breakdowns['hcap'][0]
                            away_cover = breakdowns['hcap'][2]
 
                            kelly = calc_kelly(home_cover_odds, home_cover)
                            
                            bets.append((home, kelly, hcap, home_cover_odds)) 
                            kelly = calc_kelly(away_cover_odds, away_cover)
                            bets.append((away, kelly, -hcap, away_cover_odds))
                            if game[0] == "EC":
                                orig_home_cover_odds = orig_away_cover_odds = 1.01
                            #hcap adjust
                            
                            imp = breakdowns['hcap_2']
                            if imp[1] == 0 and imp[2] == 0:
                                kelly = kelly_hl(orig_home_cover_odds, imp[0], imp[3])       
                            elif imp[2] == 0 and imp[3] == 0:
                                kelly = kelly_hw(orig_home_cover_odds, imp[0], imp[1])
                            else:
                                kelly = calc_kelly(orig_home_cover_odds, imp[0], imp[2])
                            bets.append((home, kelly, orig, orig_home_cover_odds))
                  
                            if imp[1] == 0 and imp[2] == 0:
                                kelly = kelly_hw(orig_away_cover_odds, imp[4], imp[3])
                            elif imp[2] == 0 and imp[3] == 0:
                                kelly = kelly_hl(orig_away_cover_odds, imp[4], imp[1])
                            else:
                                kelly = calc_kelly(orig_away_cover_odds, imp[4], imp[2])
                            bets.append((away, kelly, -orig, orig_away_cover_odds))
                            
                            
                            
                            #TOTALS
                            o_odds = 1 / (other_bdowns['game_total'][0] / (1 - other_bdowns['game_total'][1]) * 1.08)
                            u_odds = 1 / (other_bdowns['game_total'][2] / (1 - other_bdowns['game_total'][1]) * 1.08)
                            #odds adjust

                            game_over = breakdowns['game_total'][0]
                            game_under = breakdowns['game_total'][2]
                            kelly = calc_kelly(o_odds, game_over)
                            bets.append(('O-'+home+' vs '+away, kelly, rounded_line, o_odds))
                            kelly = calc_kelly(u_odds, game_under)
                            bets.append(('U-'+home+' vs '+away, kelly, rounded_line, u_odds))
                            total = other_bdowns['game_total_2']
                            
                            o_odds = 1 / (other_bdowns['game_total'][0] / (1 - other_bdowns['game_total'][1]) * 1.05)
                            u_odds = 1 / (other_bdowns['game_total'][2] / (1 - other_bdowns['game_total'][1]) * 1.05)
        
                            if rounded_line_2 == rounded_line:
                                o_odds_2 = o_odds
                                u_odds_2 = u_odds
                            elif rounded_line_2.is_integer():
                                o_odds_2 = 1 / (total[0] / (1 - total[2]) * 1.05)
                                u_odds_2 = (-o_odds_2) / (-o_odds_2 * 1.05 + 1)  
                            else:
                                o_whole = 1 / (total[0] / (1 - max(total[1:4])) * 1.05)
                                u_whole = (-o_whole) / (-o_whole * 1.05 + 1) 
                                o_odds_2 = (o_odds+o_whole)/2
                                u_odds_2 = (u_odds+u_whole)/2
                            #odds2 adjust 
                            
                            total = breakdowns['game_total_2']
                            if rounded_line == rounded_line_2:
                                kelly = calc_kelly(o_odds, total[0])
                                bets.append(('O-'+home+' vs '+away, kelly, rounded_line_2, o_odds_2))
                                kelly = calc_kelly(u_odds, total[4])
                                bets.append(('U-'+home+' vs '+away, kelly, rounded_line_2, u_odds_2))
                            elif rounded_line_2.is_integer():
                                kelly = calc_kelly(o_odds_2, total[0], total[2])
                                bets.append(('O-'+home+' vs '+away, kelly, rounded_line_2, o_odds_2))
                                kelly = calc_kelly(u_odds_2, total[4], total[2])
                                bets.append(('U-'+home+' vs '+away, kelly, rounded_line_2, u_odds_2))
                            elif total[1] > 0:
                                kelly = kelly_hw(o_odds_2, total[0], total[1])
                                bets.append(('O-'+home+' vs '+away, kelly, rounded_line_2, o_odds_2))
                                kelly = kelly_hl(u_odds_2, total[4], total[1])
                                bets.append(('U-'+home+' vs '+away, kelly, rounded_line_2, u_odds_2))
                            elif total[3] > 0:
                                kelly = kelly_hl(o_odds_2, total[0], total[3])
                                bets.append(('O-'+home+' vs '+away, kelly, rounded_line_2, o_odds_2))
                                kelly = kelly_hw(u_odds_2, total[4], total[3])
                                bets.append(('U-'+home+' vs '+away, kelly, rounded_line_2, u_odds_2))
                                    
                            ho_odds = 1 / (other_bdowns['home_total'][0] / (1 - other_bdowns['home_total'][1]) * 1.07)
                            hu_odds = 1 / (other_bdowns['home_total'][2] / (1 - other_bdowns['home_total'][1]) * 1.07)
                            ao_odds = 1 / (other_bdowns['away_total'][0] / (1 - other_bdowns['away_total'][1]) * 1.07)
                            au_odds = 1 / (other_bdowns['away_total'][2] / (1 - other_bdowns['away_total'][1]) * 1.07)
                            #adjust
    

                            home_over = breakdowns['home_total'][0]
                            home_under = breakdowns['home_total'][2]
                            kelly = calc_kelly(ho_odds, home_over, breakdowns['home_total'][1])
                            bets.append(('O-'+home, kelly, hrounded_line, ho_odds))
                            kelly = calc_kelly(hu_odds, home_under, breakdowns['home_total'][1])
                            bets.append(('U-'+home, kelly, hrounded_line, hu_odds))
                                   
                            away_over = breakdowns['away_total'][0]
                            away_under = breakdowns['away_total'][2]
                            kelly = calc_kelly(ao_odds, away_over, breakdowns['away_total'][1])
                            bets.append(('O-'+away, kelly, arounded_line, ao_odds)) 
                            kelly = calc_kelly(au_odds, away_under, breakdowns['away_total'][1])
                            bets.append(('U-'+away, kelly, arounded_line, au_odds)) 
                                  
                            #BTTS/0.5
                            btts_yes_odds = 1 / other_bdowns['btts'][0] / 1.07
                            btts_no_odds = 1 / other_bdowns['btts'][1] / 1.07
                            btts_yes = breakdowns['btts'][0]
                            btts_no = breakdowns['btts'][1]
                            kelly = calc_kelly(btts_yes_odds, btts_yes)
                            bets.append(('BTTS-Yes-'+home+' vs '+away, kelly, rounded_line, btts_yes_odds))
                            kelly = calc_kelly(btts_no_odds, btts_no)
                            bets.append(('BTTS-No-'+home+' vs '+away, kelly, rounded_line, btts_no_odds))
                            
                            hobttsy_odds = 1 / other_bdowns['home_or_btts'][0] / 1.1   
                            dobttsy_odds = 1 / other_bdowns['draw_or_btts'][0] / 1.1 
                            aobttsy_odds = 1 / other_bdowns['away_or_btts'][0] / 1.1 
                            hobttsn_odds = 1 / other_bdowns['home_or_btts'][1] / 1.1                      
                            dobttsn_odds = 1 / other_bdowns['draw_or_btts'][1] / 1.1 
                            aobttsn_odds = 1 / other_bdowns['away_or_btts'][1] / 1.1
                            if hrounded_line == .5:
                                hobttsy_odds = ho_odds
                            if arounded_line == .5:
                                aobttsy_odds = ao_odds
                            
                                                      
                          
                            kelly = calc_kelly(hobttsy_odds, breakdowns['home_or_btts'][0])
                            if hobttsy_odds > 1:
                                bets.append(('HOBTTSY- '+home, kelly, hcap, hobttsy_odds))
                            kelly = calc_kelly(hobttsn_odds, breakdowns['home_or_btts'][1])
                            bets.append(('HOBTTSN- '+home, kelly, hcap, hobttsn_odds))
                            kelly = calc_kelly(dobttsy_odds, breakdowns['draw_or_btts'][0])
                            bets.append(('DOBTTSY- '+home+' vs '+away, kelly, hcap, dobttsy_odds))
                            kelly = calc_kelly(dobttsn_odds, breakdowns['draw_or_btts'][1])
                            bets.append(('DOBTTSN- '+home+' vs '+away, kelly, hcap, dobttsn_odds))
                            kelly = calc_kelly(aobttsy_odds, breakdowns['away_or_btts'][0])
                            bets.append(('AOBTTSY- '+away, kelly, hcap, aobttsy_odds))
                            kelly = calc_kelly(aobttsn_odds, breakdowns['away_or_btts'][1])
                            bets.append(('AOBTTSN- '+away, kelly, hcap, aobttsn_odds))
                            
                            
                            #DNB/ANB/HNB/DC
                            hdnb_odds = 1 / (other_bdowns['result'][0] / (1 - other_bdowns['result'][1]) * 1.08) 
                            adnb_odds = 1 / (other_bdowns['result'][2] / (1 - other_bdowns['result'][1]) * 1.08)                             
                            hanb_odds = 1 / (other_bdowns['result'][0] / (1 - other_bdowns['result'][2]) * 1.1)
                            danb_odds = 1 / (other_bdowns['result'][1] / (1 - other_bdowns['result'][2]) * 1.1)
                            ahnb_odds = 1 / (other_bdowns['result'][2] / (1 - other_bdowns['result'][0]) * 1.1)
                            dhnb_odds = 1 / (other_bdowns['result'][1] / (1 - other_bdowns['result'][0]) * 1.1)
                            
                            #dnb adjust
                      


                            kelly = (home_win * (hdnb_odds-1) - (away_win))/((home_win+away_win)*(hdnb_odds-1))
                            if hdnb_odds > 1:
                                bets.append(('DNB- '+home, kelly, hcap, hdnb_odds))
                            kelly = (away_win * (adnb_odds-1) - (home_win))/((home_win+away_win)*(adnb_odds-1))
                            if adnb_odds > 1:
                                bets.append(('DNB- '+away, kelly, hcap, adnb_odds))
                                
                            kelly = (home_win * (hanb_odds-1) - (draw))/((home_win+draw)*(hanb_odds-1))
                            if hanb_odds > 1:
                                bets.append(('ANB- '+home, kelly, hcap, hanb_odds))
                
                            kelly = (draw * (danb_odds-1) - (home_win))/((home_win+draw)*(danb_odds-1))
                            if danb_odds > 1:
                                bets.append(('ANB- DRAW-'+home+' '+ away, kelly, hcap,danb_odds))
                               
                            kelly = (away_win * (ahnb_odds-1) - (draw))/((away_win+draw)*(ahnb_odds-1))
                            if ahnb_odds > 1:
                                bets.append(('HNB- '+away, kelly, hcap, ahnb_odds))
                
                            kelly = (draw * (dhnb_odds-1) - (away_win))/((away_win+draw)*(dhnb_odds-1))
                            if dhnb_odds > 1:
                                bets.append(('HNB- DRAW-'+home+' '+ away, kelly, hcap,dhnb_odds))
                                
                            hdc_odds = 1 / (other_bdowns['result'][0] + other_bdowns['result'][1]) / 1.1
                            adc_odds = 1 / (other_bdowns['result'][2] + other_bdowns['result'][1]) / 1.1
                            hdc = breakdowns['result'][0] + breakdowns['result'][1]
                            adc = breakdowns['result'][2] + breakdowns['result'][1]
                            kelly = calc_kelly(hdc_odds, hdc)
                            if hdc_odds > 1:
                                bets.append(('DC- '+home, kelly, hcap, hdc_odds))
                            kelly = calc_kelly(adc_odds, adc)
                            if adc_odds > 1:
                                bets.append(('DC- '+away, kelly, hcap, adc_odds))
                                                
                            #RESULT+BTTS
                            hbttsy_odds = 1 / other_bdowns['win+btts'][0] / 1.4
                            abttsy_odds = 1 / other_bdowns['win+btts'][1] / 1.4 
                            dbttsy_odds = 1 / 1 / other_bdowns['draw+btts'][0] / 1.4 
                            dbttsn_odds = 1 / other_bdowns['draw+btts'][1] / 1.4
                            kelly = calc_kelly(hbttsy_odds, breakdowns['win+btts'][0])
                            bets.append(('HBTTSY- '+home, kelly, hcap, hbttsy_odds))
                            kelly = calc_kelly(abttsy_odds, breakdowns['win+btts'][1])
                            bets.append(('ABTTSY- '+away, kelly, hcap, abttsy_odds))
                            kelly = calc_kelly(dbttsy_odds, breakdowns['draw+btts'][0])
                            bets.append(('DBTTSY- '+home+' vs '+away, kelly, hcap, dbttsy_odds))
                            kelly = calc_kelly(dbttsn_odds, breakdowns['draw+btts'][1])
                            bets.append(('DBTTSN- '+home+' vs '+away, kelly, hcap, dbttsn_odds))
                            htny_odds = 1 / other_bdowns['home_to_nil'][0] / 1.1
                            htnn_odds = 1 / other_bdowns['home_to_nil'][1]  / 1.1                 
                            atny_odds = 1 / other_bdowns['away_to_nil'][0] / 1.1
                            atnn_odds = 1 / other_bdowns['away_to_nil'][1] / 1.1
                            if atnn_odds <= 1:
                                atnn_odds = 1.01

                         
                            kelly = calc_kelly(htny_odds, breakdowns['home_to_nil'][0])
                            bets.append(('HTNY- '+home+' vs '+away, kelly, hcap, htny_odds))
                            kelly = calc_kelly(htnn_odds, breakdowns['home_to_nil'][1])
                            bets.append(('HTNN- '+home+' vs '+away, kelly, hcap, htnn_odds))
                            kelly = calc_kelly(atny_odds, breakdowns['away_to_nil'][0])
                            bets.append(('ATNY- '+home+' vs '+away, kelly, hcap, atny_odds))
                            kelly = calc_kelly(atnn_odds, breakdowns['away_to_nil'][1])
                            bets.append(('ATNN- '+home+' vs '+away, kelly, hcap, atnn_odds))
                                                    
                            best = best2 = best3 = -float('inf')
                            pick = pick2 = pick3 = None
                            
                            
                            for bet in bets:
                                if let == 0:
                                    betlists[home+away] = betlists.setdefault(home+away, []) + [bet]        
                                if bet[1] > best and bet[3]>=1.5 and bet[3]<=3:
                                    best = bet[1]
                                    pick = bet
                                if bet[1] > best2 and bet[3]>=3 and bet[3]<=5:
                                    best2 = bet[1]
                                    pick2 = bet
                                if bet[1] > best3 and bet[3]>=down and bet[3]<=up:            
                                    best3 = bet[1]
                                    pick3 = bet
                            #sort(betlists[home+away], 0, .03759)
                                    
                            if pick and pick[1] >= -1 and pick[1] <= 1 and let == 0:
                                
                                one.append(other_bdowns)
                                two.append(breakdowns)
                                picks.append(pick)
                                check.append(home)
                                #for bet in bets:
                                #    if bet[1] >= .02929 and bet[1] <= .0421 and let == 0 and bet[3]>=1.5 and bet[3]<=3:
                                #            print(bet)
                            #.00954 .01123                   
                            if pick2 and pick2[1] >= -1 and pick2[1] <= 1 and let == 1:
                                dog_picks.append(pick2)
                                
                            #.1081, .1272
                            if pick3 and pick3[1] >= -1 and pick3[1] <= 1 and let == 2:
                                fav_picks.append(pick3)
                                
                            if len(div_memo) >= 2000000:
                                div_memo = {}
                                
                            if len(mult_memo) >= 2000000:
                                mult_memo = {}
                                
                            if len(prob_memo) >= 2000000:
                                prob_memo = {}
                                
                            if len(exp_memo) >= 2000000:
                                exp_memo = {}
                        

        


#sort(fav_picks, 0, .12218)
#sort(dog_picks, 0, .00979)
#sort(picks, 0, .03759)

#print(draw_picks)
#print(picks)
#print(dog_picks)
#print(len(picks))
#print(check)


