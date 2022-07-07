from time import sleep


def get_tree_fpbase(driver, root):
    driver.get("https://www.fpbase.org/protein/" + root)
    sleep(1)

    return driver.execute_script(
        "return [...document.querySelectorAll('.lineage-wrapper svg')[0].getElementsByClassName('node')].map(v=>["
        "v.textContent, v.getElementsByTagName('a')[0].href.baseVal])")


def crawl_dataset(tar):
    import pandas as pd
    return pd.read_html("https://www.fpbase.org" + tar)


def crawl_data(driver, tar):
    from selenium.webdriver.common.by import By
    driver.get("https://www.fpbase.org" + tar)
    sleep(0.5)

    return driver.execute_script("return document.querySelector('.formatted_aminosquence').textContent."
                                 "replaceAll(' ','')"), driver.find_element(By.CLASS_NAME, "reference-details").find_element(By.CSS_SELECTOR, "a").get_attribute("href")
