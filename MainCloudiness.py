from Weather.Cloudiness import CloudCalculator

if __name__ == '__main__':
    area_center_lat = 55.365
    area_center_long = 38.305
    area_range_lat = 2
    area_range_long = 2
    cloud_calculator = CloudCalculator()
    cloud_calculator.to_set_coordinates_of_calculation_area(area_center_lat, area_center_long,
                                                            area_range_lat, area_range_long)
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2007.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2008.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2009.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2010.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2011.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2012.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2013.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2014.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2015.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2016.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2017.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2018.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2019.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2020.nc')
    cloud_calculator.to_add_cloud_dataset('D:/Data/NetCDF cloud data/Total cloud cover/tcdc.eatm.gauss.2021.nc')

    cloud_distr_dataset = cloud_calculator.to_calculate_cloud_distr()
    cloud_distr_dataset.to_save_data('Бронницкое_лесничество_cloud_cover', 'D:/Проекты/Оперативность/Бронницкое лесничество/Промежуточные данные')