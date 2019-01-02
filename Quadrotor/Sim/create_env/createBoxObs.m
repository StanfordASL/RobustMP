function obs = createBoxObs(origin, edges)
    obs = [origin; origin+edges];
end